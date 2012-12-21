# DEMO naive Bayes Classification
library(predictiveModeling)
library(class)
library(e1071)
library(biclust)
library(clinfun)
library(randomForest)
library(synapseClient)
library(multicore)
library(doMC)
registerDoMC()

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModelLogistic.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myNaiveBayesModel.R")


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


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

Folder = "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/"
for(k1 in 1:ncol(adf_drug)){
  
  RFList<-c()
  
  filename = paste(Folder,"randomForest/ROC_",k1,".Rdata",sep="")
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,k1,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  RFList<-crossValidatePredictiveModelLogistic(filteredFeatureDataScaled,filteredResponseDataScaled,model = myRandomForestModel$new())  
  save(RFList,file=filename)
  #NBList[[k1]]<-crossValidatePredictiveModelLogistic(filteredFeatureDataScaled,filteredResponseDataScaled,model = myNaiveBayesModel$new())  
  
  print(k1)
}
