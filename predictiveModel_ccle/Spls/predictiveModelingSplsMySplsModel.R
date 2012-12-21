### DEMO Ridge
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/mySplsModel.R")
###################################################
#### Load Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "48339"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$copySet

id_oncomapLayer <- "48341"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$oncomapSet

id_exprLayer <- "48344" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$exprSet

###################################################
### Load Response Data
###################################################

id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


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
  
  etas <- seq(0.1,0.9,0.1)
  Ks <- seq(1,10,by =1)
  
  
  resultsOriginal<-crossValidatePredictiveModel1(filteredFeatureDataScaled,filteredResponseDataScaled, model = mySplsModel$new(), eta = etas, K = Ks, nfolds = 5)
  
  filename <- paste("cvDrug_",kk,"_Spls.Rdata",sep ="")  
  save(resultsOriginal,lambdas,file = filename)
}
