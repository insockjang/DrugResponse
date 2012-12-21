### DEMO ENet
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical1.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")
# library(multicore)
# library(doMC)
# registerDoMC()

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

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

for(kk in 1:ncol(pData(adf_drug))){
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)    
  
  resultsCat<-crossValidatePredictiveModel_categorical1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myRandomForestModel$new(), numFolds = 5,ntree = 100, thresholdMethod = "median")
  filename = paste("randomForest/ROCR_median_",kk,".Rdata",sep = "")
  save(resultsCat,file = filename)
  
  resultsCat<-crossValidatePredictiveModel_categorical1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myRandomForestModel$new(), numFolds = 5,ntree = 100, thresholdMethod = "median_mad")
  filename = paste("randomForest/ROCR_median_mad_",kk,".Rdata",sep = "")
  save(resultsCat,file = filename)
  
  resultsCat<-crossValidatePredictiveModel_categorical1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myRandomForestModel$new(), numFolds = 5,ntree = 100, thresholdMethod = "mean")
  filename = paste("randomForest/ROCR_mean_",kk,".Rdata",sep = "")
  save(resultsCat,file = filename)
  
  resultsCat<-crossValidatePredictiveModel_categorical1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myRandomForestModel$new(), numFolds = 5,ntree = 100, thresholdMethod = "mean_sd")
  filename = paste("randomForest/ROCR_mean_sd_",kk,".Rdata",sep = "")
  save(resultsCat,file = filename)
  
  print(kk)  
}

