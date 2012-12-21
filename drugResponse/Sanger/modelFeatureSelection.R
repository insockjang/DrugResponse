### DEMO ENet
library(multicore)
library(doMC)
registerDoMC()

library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")


####################################### Sanger
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


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

rankTransform <- function(FEA){
  FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
}


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
  
  alphas =unique(createENetTuneGrid()[,1])
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  modelLasso <- myEnetModel$new()
  modelLasso$customTrain(filteredFeatureDataScaled, filteredResponseDataScaled, alpha=1, lambda = lambdas, nfolds = 5)
  resultsLasso<-modelLasso$getCoefficients()
  
  modelRidge <- myEnetModel$new()
  modelRidge$customTrain(filteredFeatureDataScaled, filteredResponseDataScaled, alpha=10^-8, lambda = lambdas, nfolds = 5)
  resultsRidge<-modelRidge$getCoefficients()
  
  modelENet <- myEnetModel$new()
  modelENet$customTrain(filteredFeatureDataScaled, filteredResponseDataScaled, alpha=alphas, lambda = lambdas, nfolds = 5)
  resultsENet<-modelENet$getCoefficients()
  
  filename1 <- paste("Lasso2/fsDrug_",kk,".Rdata",sep ="")  
  filename2 <- paste("Ridge2/fsDrug_",kk,".Rdata",sep ="")  
  filename3 <- paste("ENet2/fsDrug_",kk,".Rdata",sep ="")  
  save(resultsLasso,file = filename1)
  save(resultsRidge,file = filename2)
  save(resultsENet,file = filename3)
}
