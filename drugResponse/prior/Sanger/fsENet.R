
library(predictiveModeling)
library(synapseClient)
library(multicore)
library(doMC)
registerDoMC()

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/bootstrapPredictiveModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")

source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("~/COMPBIO/trunk/users/jang/drugResponse/priorFilterData.R")
source("~/COMPBIO/trunk/users/jang/drugResponse/filterPredictiveModelData1.R")


###################################################
#### Load Sanger Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "210937"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "266141"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_exprLayer <- "210931" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "220680" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

id_priorLayer <- "273739"
layer_prior <- loadEntity(id_priorLayer)
priorListCGC <- layer_prior$objects$cancerGeneCensus
priorListTF <- layer_prior$objects$transcriptionFactor
priorListM <- layer_prior$objects$modulator
priorListR <- union(layer_prior$objects$transcriptionFactor,layer_prior$objects$modulator)

dataSet1<-priorFilterData(eSet_expr, eSet_copy, eSet_oncomap,adf_drug,priorList=priorListCGC)
dataSet2<-priorFilterData(eSet_expr, eSet_copy, eSet_oncomap,adf_drug,priorList=priorListTF)
dataSet3<-priorFilterData(eSet_expr, eSet_copy, eSet_oncomap,adf_drug,priorList=priorListM)
dataSet4<-priorFilterData(eSet_expr, eSet_copy, eSet_oncomap,adf_drug,priorList=priorListR)

for(kk in 1:ncol(dataSet1$responseData)){  
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
    Data1<-filterPredictiveModelData1(dataSet1, kk, refineFilter = TRUE)
  Data2<-filterPredictiveModelData1(dataSet2, kk, refineFilter = TRUE)
  Data3<-filterPredictiveModelData1(dataSet3, kk, refineFilter = TRUE)
  Data4<-filterPredictiveModelData1(dataSet4, kk, refineFilter = TRUE)
  
  
  # 5 fold cross validation 
#   alphas <- 10^-8
  #   lambdas <- exp(seq(-5, 2, length = 100))
  alphas = unique(createENetTuneGrid()[,1])
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
#   resultsCGC<-crossValidatePredictiveModel1(Data1$featureData, Data1$responseData, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
#   resultsTF<-crossValidatePredictiveModel1(Data2$featureData, Data2$responseData, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
#   resultsM<-crossValidatePredictiveModel1(Data3$featureData, Data3$responseData, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
#   resultsR<-crossValidatePredictiveModel1(Data4$featureData, Data4$responseData, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
  
  
  enetCGC <- myEnetModel$new()
  enetCGC$customTrain(Data1$featureData, Data1$responseData, alpha = alphas, lambda = lambdas, nfolds = 5)
  resultsCGC<-enetCGC$getCoefficients()
  
  enetTF <- myEnetModel$new()
  enetTF$customTrain(Data2$featureData, Data2$responseData, alpha = alphas, lambda = lambdas, nfolds = 5)
  resultsTF<-enetTF$getCoefficients()
  
  enetM <- myEnetModel$new()
  enetM$customTrain(Data3$featureData, Data3$responseData, alpha = alphas, lambda = lambdas, nfolds = 5)
  resultsM<-enetM$getCoefficients()
  
  enetR <- myEnetModel$new()
  enetR$customTrain(Data4$featureData, Data4$responseData, alpha = alphas, lambda = lambdas, nfolds = 5)
  resultsR<-enetR$getCoefficients()
  
  filename <- paste("ENet2/fsDrug_",kk,".Rdata",sep ="")  
  
  save(resultsCGC,resultsTF,resultsM,resultsR,file = filename)  
}
