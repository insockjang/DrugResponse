### DEMO bootstrapping
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
  Data1<-filterPredictiveModelData1(dataSet1, kk, refineFilter = FALSE)
  Data2<-filterPredictiveModelData1(dataSet2, kk, refineFilter = FALSE)
  Data3<-filterPredictiveModelData1(dataSet3, kk, refineFilter = FALSE)
  Data4<-filterPredictiveModelData1(dataSet4, kk, refineFilter = FALSE)
  
  
  # 5 fold cross validation 
  alphas =1
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  resultsCGC<-bootstrapPredictiveModel(Data1$featureData, Data1$responseData, model = myEnetModel$new(), numBootstrap=100,alpha=alphas, lambda = lambdas)
  resultsTF<-bootstrapPredictiveModel(Data2$featureData, Data2$responseData, model = myEnetModel$new(), numBootstrap=100,alpha=alphas, lambda = lambdas)
  resultsM<-bootstrapPredictiveModel(Data3$featureData, Data3$responseData, model = myEnetModel$new(), numBootstrap=100,alpha=alphas, lambda = lambdas)
  resultsR<-bootstrapPredictiveModel(Data4$featureData, Data4$responseData, model = myEnetModel$new(), numBootstrap=100,alpha=alphas, lambda = lambdas)
  
  filename <- paste("bootstrapLasso2/cvDrug_",kk,".Rdata",sep ="")  
  save(resultsCGC,resultsTF,resultsM,resultsR,file = filename)  
}
