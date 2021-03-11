
library(randomForest)
library(pROC)
library(graphics)

OTU_Meta <- read.csv("~/41467_2020_14677_MOESM5_ESM_cleantoprows.csv", header = T)
OTU_Meta = OTU_Meta[1:264,]


####### Sort the features based on their averaged mean decreased Gini score, which is averaged from total 264 leave-one-out test.
#make global data frame for importance scores (gini scores)
OTUginiAccuracy = data.frame("taxa" = colnames(OTU_Meta)[-(1:3)])
##run model on all training sets (leave one out)
## takes validation sets using each of 264 patients (1 patient per set "leave one out"), iterate through all 264
## training sets composed of 263 patients 
for (i in 1:nrow(OTU_Meta[,-1])) {
  
  ValidSet = OTU_Meta[i,-(1:2)]
  TrainSet = OTU_Meta[-i,-(1:2)]
  # factorize the response variable
  TrainSet$CST <- factor(TrainSet$CST)
  ValidSet$CST <- factor(ValidSet$CST)
  
  # Create a Random Forest model with default parameters
  OTU_Model1 <- randomForest(CST ~ ., data = TrainSet, importance = TRUE)
  
  #finds the gini scores using the importance function
  OTUginiAccuracyTemp = as.data.frame(importance(OTU_Model1, type=2))
  
  ##record the gini scores into the global data frame above
  OTUginiAccuracy = cbind(OTUginiAccuracy, OTUginiAccuracyTemp$MeanDecreaseGini)
}

#average gini scores of each taxa (collected from 264 iterations of training sets)
meanAccuracy_OTU=rowMeans(OTUginiAccuracy[,-1], na.rm = FALSE, dims = 1)
meanAccuracy_OTUtable = data.frame(OTUginiAccuracy$taxa, as.numeric(as.character(meanAccuracy_OTU)))
colnames(meanAccuracy_OTUtable) =c("Taxa", "mean_decreaseGigi")
#sort gini scores and matching taxa names, obtain list of ranked taxa based on importance
meanAccuracy_OTUtableSort = meanAccuracy_OTUtable[order(-meanAccuracy_OTUtable$mean_decreaseGigi),]
#### feature list done


######### calculate accuracy of the prediction using taxa sets from the ranked list
###each round, carry out leave-one-out test
##"rank" = number of top taxa used (e.g. rank = 1 uses 1 top taxa; rank = 2 uses 2 top taxa; rank = 28 uses 28 top taxa)
EachRankAccuracy_Results = list()
RankAccuracyRatios = list()

#iteration; recursively adds the next taxa to the "set" in each rank, starting from the first (most important) taxa
### finds the PREDICTION ACCURACY based on VALIDATION SETS
##parse taxa name -- allows it to stay as a data frame even when only 1 starting taxa

#To save time, I did not use all taxa, but the top 60 taxa as the upper limit
for (j in 2:60) {
  Feature = c("CST", as.vector(meanAccuracy_OTUtableSort$Taxa[1:j]))
  
  #add metadata back (response status)
  NewOTU_Meta = as.data.frame(OTU_Meta[,colnames(OTU_Meta) %in% Feature])
  print (j)
  summary = list()
  
  # use cross validation to find the prediction accuracy
  for (i in 1:nrow(NewOTU_Meta)) {
    print (i)
    ### same leave one out sets as before
    ValidSet = NewOTU_Meta[i,]
    TrainSet = NewOTU_Meta[-i,]
    #reassign factors
    TrainSet$CST <- factor(TrainSet$CST)
    ValidSet$CST <- factor(ValidSet$CST)
    
    # Create a Random Forest model with default parameters
    NewOTU_Model <- randomForest(CST ~ ., data = TrainSet, importance = TRUE)
    
    predValid <- predict(NewOTU_Model, ValidSet, type = "prob")
    predValid_class <- predict(NewOTU_Model, ValidSet, type = "class")

    ## Checking classification accuracy
    if(as.character(predValid_class) == as.character(ValidSet[,1])){
      isCorrect = 1
    }else{
      isCorrect = 0
    }
    
    #### Prediction rate
    summary[[i]] = as.matrix(t(c(i, j, NewOTU_Model$err.rate[nrow(NewOTU_Model$err.rate),1],predValid[1],predValid[2], predValid[3], predValid[4],predValid[5], predValid[6], as.character(ValidSet[,1]), isCorrect)))
    
  } #end of leave-one-out
  EachRankAccuracy_Results[[j]] = data.frame(t(sapply(summary, "[")))
  
  ##sum of isCorrect (correct prediction total)
  correctTotal = sum(as.numeric(as.character(EachRankAccuracy_Results[[j]][,7])))
  #find fraction out of total patients possible
  correctRatio = correctTotal/41
  
  #error average
  errorAverage = mean(as.numeric(as.character(EachRankAccuracy_Results[[j]][,3])))
  RankAccuracyRatios[[j]] = as.matrix(t(c(j,correctRatio, errorAverage)))
  
}
###### a summary dataframe!
### showcases: sample#, rank (# of taxa from top), OOB error rate, CST I prediction rate, CST II prediction rate, CST III prediction rate,CST IV-A prediction rate,CST IV-B prediction rate, CST V prediction rate,nonresponder label (if guessed nonresponder), the true status, and whether the prediction was correct (1 or 0)
totalsummary = data.frame(do.call(rbind,EachRankAccuracy_Results))
colnames(totalsummary) = c("sample", "rank", "OOB", "I","II","III", "IV-A","IV-B","V", "Truth", "Correctness")


######### Calculate AUC

HighestRank_allRow <- max(as.numeric(totalsummary$rank))
auc_res <- list()

#iteration for each set/rank of taxa -- finds AUC for each iteration
for (NumberOfFeature in 1:HighestRank_allRow) {
  tryCatch({
    RankData <- totalsummary[totalsummary$rank == NumberOfFeature,]
    predictionrate=data.frame(RankData$I,RankData$II,RankData$III, RankData$`IV-A`, RankData$`IV-B`, RankData$V)
    colnames(predictionrate) = c("I", "II", "III", "IV-A", "IV-B", "V")
    AUC_Area=auc(multiclass.roc(RankData$Truth, predictionrate))
    auc_res[[NumberOfFeature]] <- data.frame(RankData[1,1:2], NumberOfFeature, AUC_Area)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

####data frame of auc
##features: sample#, rank, NumberOfFeature (same thing as rank), AUC Area
auc = data.frame(do.call(rbind,auc_res))



############################# ROC curve using Rank= 4
###Pick the best rank = 4
TopRankSummary = totalsummary[totalsummary$rank ==4,]
predictionrate=data.frame(as.numeric(as.character(TopRankSummary$I)),as.numeric(as.character(TopRankSummary$II)),as.numeric(as.character(TopRankSummary$III)), as.numeric(as.character(TopRankSummary$`IV-A`)), as.numeric(as.character(TopRankSummary$`IV-B`)), as.numeric(as.character(TopRankSummary$V)))
colnames(predictionrate) = c("I", "II", "III", "IV-A", "IV-B", "V")
AUC_Area=auc(multiclass.roc(TopRankSummary$Truth, predictionrate))


for (i in 1:nrow(TopRankSummary)) {
  if(TopRankSummary$Truth[i] == "I") {
    TopRankSummary$I_label[i] = 1
    TopRankSummary$II_label[i] = 0
    TopRankSummary$III_label[i] = 0
    TopRankSummary$IV_A_label[i] = 0
    TopRankSummary$IV_B_label[i] = 0
    TopRankSummary$V_label[i] = 0
  }else if(TopRankSummary$Truth[i] == "II") {
    TopRankSummary$I_label[i] = 0
    TopRankSummary$II_label[i] = 1
    TopRankSummary$III_label[i] = 0
    TopRankSummary$IV_A_label[i] = 0
    TopRankSummary$IV_B_label[i] = 0
    TopRankSummary$V_label[i] = 0
  }else if(TopRankSummary$Truth[i] == "III") {
    TopRankSummary$I_label[i] = 0
    TopRankSummary$II_label[i] = 0
    TopRankSummary$III_label[i] = 1
    TopRankSummary$IV_A_label[i] = 0
    TopRankSummary$IV_B_label[i] = 0
    TopRankSummary$V_label[i] = 0
  }else if(TopRankSummary$Truth[i] == "IV-A") {
    TopRankSummary$I_label[i] = 0
    TopRankSummary$II_label[i] = 0
    TopRankSummary$III_label[i] = 0
    TopRankSummary$IV_A_label[i] = 1
    TopRankSummary$IV_B_label[i] = 0
    TopRankSummary$V_label[i] = 0
  }else if(TopRankSummary$Truth[i] == "IV-B") {
    TopRankSummary$I_label[i] = 0
    TopRankSummary$II_label[i] = 0
    TopRankSummary$III_label[i] = 0
    TopRankSummary$IV_A_label[i] = 0
    TopRankSummary$IV_B_label[i] = 1
    TopRankSummary$V_label[i] = 0
  }else {
    TopRankSummary$I_label[i] = 0
    TopRankSummary$II_label[i] = 0
    TopRankSummary$III_label[i] = 0
    TopRankSummary$IV_A_label[i] = 0
    TopRankSummary$IV_B_label[i] = 0
    TopRankSummary$V_label[i] = 1
  }
}

ROC1= roc(as.factor(TopRankSummary$I_label), as.numeric(as.vector(TopRankSummary$I)))
ROC2= roc(as.factor(TopRankSummary$II_label), as.numeric(as.vector(TopRankSummary$II)))
ROC3= roc(as.factor(TopRankSummary$III_label), as.numeric(as.vector(TopRankSummary$III)))
ROC4A= roc(as.factor(TopRankSummary$IV_A_label), as.numeric(as.vector(TopRankSummary$`IV-A`)))
ROC4B= roc(as.factor(TopRankSummary$IV_B_label), as.numeric(as.vector(TopRankSummary$`IV-B`)))
ROC5= roc(as.factor(TopRankSummary$V_label), as.numeric(as.vector(TopRankSummary$V)))

plot(ROC1, col = 1, lty = 3, main = "ROC")  #black
plot(ROC2, col = 2, lty = 3, add = TRUE)  #red
plot(ROC3, col = 3, lty = 3, add = TRUE)   #green
plot(ROC4A, col = 4, lty = 3, add = TRUE)   #blue
plot(ROC4B, col = 5, lty = 3, add = TRUE)  #light blue
plot(ROC5, col = 6, lty = 3, add = TRUE)  #purple
