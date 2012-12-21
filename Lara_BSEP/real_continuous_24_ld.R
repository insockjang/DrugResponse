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
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")

###################################################
#### Load BSEP Molecular Feature Data from Synapse ####
###################################################

data_expr <- train24ld
data_drug<-log(as.matrix(train24IC50))


# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])
filteredFeatureData  <- (filteredData$featureData)
filteredResponseData <- (filteredData$responseData)

alphas = seq(0.1,1,by = 0.1)
lambdas = createENetTuneGrid(alphas = 1)[,2]


modelLasso<-myEnetModel$new()
modelLasso$customTrain(filteredFeatureData,filteredResponseData,alpha=1,lambda = lambdas, nfolds=5)

modelENet<-myEnetModel$new()
modelENet$customTrain(filteredFeatureData,filteredResponseData,alpha=alphas,lambda = lambdas, nfolds=5)

modelRFCat<-myRandomForestModel$new()
modelRFCat$customTrain(filteredFeatureData,filteredResponseData,ntrees = 300)

fit<-randomForest(filteredFeatureData,filteredResponseData,ntrees = 300)

AA<-as.matrix(test24ld[,which(is.na(apply(test24fc,2,mean))==0)])
predicted<-predict(fit,t(AA),type = "response")
observed <- as.numeric(test24IC50[which(is.na(apply(test24fc,2,mean))==0)])


allTePred_Lasso <- modelLasso$customPredict((t(as.matrix(test24ld))))
allTePred_ENet <- modelENet$customPredict((t(as.matrix(test24ld))))
allTePred_RF <- predicted

allTeObsr <- as.numeric(test24IC50)


catenateResults1<-function(allTePred,allTeObsr){
  TestPr<-allTePred
  TestOb<-allTeObsr
  return(list(cor = cor(TestPr,TestOb,use="complete"),
              predicted = TestPr,
              observed = TestOb))
}

resultsLasso<-catenateResults1(allTePred_Lasso,allTeObsr)
resultsENet<-catenateResults1(allTePred_ENet,allTeObsr)
resultsRandomForest_300<-catenateResults1(allTePred_RF,allTeObsr)

filename = paste("~/COMPBIO/trunk/users/jang/Lara_BSEP/ROCR_continuous_real_24ld.Rdata",sep = "")
save(resultsLasso,resultsENet,resultsRandomForest_300,file = filename)



