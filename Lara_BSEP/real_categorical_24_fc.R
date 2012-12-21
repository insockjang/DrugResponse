### LARA BSEP classifier 
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
require(multicore)
require(doMC)
registerDoMC()
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical.R")
source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")

###################################################
#### Load BSEP Molecular Feature Data from Synapse ####
###################################################

data_expr <- train24fc
Name<-names(train24catBSEP)
data_drug<-as.numeric(as.matrix(train24catBSEP))
names(data_drug)<-Name


# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])
filteredFeatureData  <- (filteredData$featureData)
filteredResponseData <- (filteredData$responseData)

alphas = seq(0.1,1,by = 0.1)
lambdas = createENetTuneGrid(alphas = 1)[,2]


modelLassoCat<-myCatEnetModel$new()
modelLassoCat$customTrain(filteredFeatureData,filteredResponseData,alpha=1,lambda = lambdas, nfolds=5)

modelENetCat<-myCatEnetModel$new()
modelENetCat$customTrain(filteredFeatureData,filteredResponseData,alpha=alphas,lambda = lambdas, nfolds=5)

modelRFCat<-myRandomForestModel$new()
modelRFCat$customTrain(filteredFeatureData,filteredResponseData,ntrees = 300)

fit<-randomForest(filteredFeatureData,filteredResponseData,ntrees = 300)


AA<-as.matrix(test24ld[,which(is.na(apply(test24fc,2,mean))==0)])
predicted<-predict(fit,t(AA),type = "response")
observed <- as.numeric(test24catBSEP[which(is.na(apply(test24fc,2,mean))==0)])


allTePred_Lasso <- modelLassoCat$customPredict((t(as.matrix(test24ld))))
allTePred_ENet <- modelENetCat$customPredict((t(as.matrix(test24ld))))
allTePred_RF <- modelRFCat$customPredict(t(test24ld))

allTeObsr <- as.numeric(test24catBSEP)



ROCresults<-function(AllTePred,AllTeObsr){
  
  a1<-which(!is.na(AllTePred)==1)  
  b<-which(!is.na(AllTeObsr)==1)  
  A1<-intersect(a1,b)  
  AllTePred<-AllTePred[A1]  
  AllTeObsr<-AllTeObsr[A1]
  
  # EVALUATE VALIDATION MODEL PERFORMANCE
  erPred <- prediction(AllTePred,AllTeObsr-1) # factor become c(1,2) from c(0,1) after concatenate
  erPerf <- performance(erPred, "tpr", "fpr")
  erAUC <- performance(erPred, "auc")
  
  # FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
  erRFPerf <- performance(erPred, "sens", "spec")
  youdensJ <- erRFPerf@x.values[[1]] + erRFPerf@y.values[[1]] - 1
  jMax <- which.max(youdensJ)
  optCut <- erPerf@alpha.values[[1]][jMax]
  
  optSens <- unlist(erRFPerf@x.values)[jMax]
  optSpec <- unlist(erRFPerf@y.values)[jMax]
  
  #     rankSum <- wilcox.test(validScoreHat[validScore == 0],validScoreHat[validScore == 1])
  
  dfPerf <- as.data.frame(cbind(unlist(erPerf@x.values), unlist(erPerf@y.values)))
  colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")
  
  rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
    geom_line() + 
    geom_abline(slope = 1, colour = "red") +
    opts(title = "Cross Validation ROC Curve") +
    ylab("True Positive Rate") +
    xlab("False Positive Rate") +
    opts(plot.title = theme_text(size = 14))
  
  return(list("erPred"=erPred,
              "erPerf"=erPerf,
              "erAUC"=erAUC,
              "rocCurve" = rocCurve))
}

resultsLasso<-ROCresults(allTePred_Lasso,allTeObsr)
resultsENet<-ROCresults(allTePred_ENet,allTeObsr)
resultsRandomForest<-ROCresults(predicted,observed)


filename = paste("~/COMPBIO/trunk/users/jang/Lara_BSEP/ROCR_categorical_real_24_fc.Rdata",sep = "")
save(resultsLasso,resultsENet,resultsRandomForest,file = filename)

