### DEMO ENet
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical1.R")

source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")
# library(multicore)
# library(doMC)
# registerDoMC()
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




featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy))

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
  
  alphas =unique(createENetTuneGrid()[,1])
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  resultsCat<-crossValidatePredictiveModel_categorical1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myCatEnetModel$new(), alpha=1, lambda = lambdas, thresholdMethod = "mean_sd")
  
  filename = paste("~/COMPBIO/trunk/users/jang/drugResponse_woMut/CCLE/catLasso2/ROCR_mean_sd_",kk,".Rdata",sep = "")
  save(resultsCat,file = filename)
  print(kk)  
}
