### DEMO bootstrapping
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/bootstrapPredictiveModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
###################################################
#### Load CCLE Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "269019"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "269021"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_exprLayer <- "269056" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "269024" 
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
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE])
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  
  # 5 fold cross validation 
  alphas =unique(createENetTuneGrid()[,1])
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  
  resultsENet<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)
  #resultsLasso<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=1, lambda = lambdas)
  #resultsRidge<-bootstrapPredictiveModel(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), numBootstrap=100, alpha=10^-5, lambda = lambdas)
  filename <- paste("bootstrapENet2/cvDrug_",kk,".Rdata",sep ="")  
  save(resultsENet,file = filename)
  
}

# ResultLasso<-c()
# ResultRidge<-c()
# ResultENet<-c()
# for(k in 1:length(resultsLasso)){
#   ResultLasso<-cbind(ResultLasso,as.matrix(resultsLasso[[k]]))
#   ResultENet<-cbind(ResultENet,as.matrix(resultsENet[[k]]))
#   ResultRidge<-cbind(ResultRidge,as.matrix(resultsRidge[[k]]))
# }
