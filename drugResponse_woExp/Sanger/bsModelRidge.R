### DEMO bootstrapping
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/bootstrapPredictiveModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")

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



featureData <- createAggregateFeatureDataSet(list(copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

Results<-list()

for(kk in 1:ncol(dataSets_ccle$responseData)){  
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  
  # 5 fold cross validation 
  alphas = 10^-8
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  
  #resultsENet<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)
  #resultsLasso<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=1, lambda = lambdas)
  resultsRidge<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)
  filename <- paste("bootstrapRidge2/cvDrug_",kk,".Rdata",sep ="")  
  save(resultsRidge,file = filename)
  
}

# ResultLasso<-c()
# ResultRidge<-c()
# ResultENet<-c()
# for(k in 1:length(resultsLasso)){
#   ResultLasso<-cbind(ResultLasso,as.matrix(resultsLasso[[k]]))
#   ResultENet<-cbind(ResultENet,as.matrix(resultsENet[[k]]))
#   ResultRidge<-cbind(ResultRidge,as.matrix(resultsRidge[[k]]))
# }
