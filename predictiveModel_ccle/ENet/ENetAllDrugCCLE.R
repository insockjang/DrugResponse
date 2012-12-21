### DEMO ENetModel + CrossValidationParameterOptimizer (alpha_vec, lambda_vec)
library(predictiveModeling)
library(synapseClient)
#library(survival)
### need to source in this file since it is not in the predictiveModeling package
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/crossValidationParameterOptimizer_penalty.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/ENetModel.R")

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


grids<-createENetTuneGrid()
alpha <- seq(0.1,1,0.1)
lambda <- createENetTuneGrid(alphas=1)[,2]

#########################################################################################################
########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
#########################################################################################################

AllDrugPredictModelsNew<-foreach(kk = 1:ncol(dataSets_ccle$responseData)) %dopar% {
  
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  featureData_processed <- filteredData$featureData
  responseData_processed <- filteredData$responseData
  
  ## scale this data
  featureData_scaled <- scale(featureData_processed[,1:500])
  responseData_scaled <- scale(responseData_processed)
  
  cv<-crossValidationParameterOptimizer_penalty$new()
  cv$findOpt(model = ENetModel,featureData_scaled,responseData_scaled,numFolds = 5,alpha = alpha, lambda = lambda)
  cv1<-crossValidationParameterOptimizer_penalty$new()
  cv1$findOpt(model = ENetModel,featureData_processed[,1:500],responseData_processed,numFolds = 10,alpha = alpha, lambda = lambda)
  
  return(list(scale=cv,nonscale=cv1))
}

#########################################################################################################
########  Analysis Step: Training and Testing data are scaled(normalized) vs. raw(unnormalized)  ########
#########################################################################################################


analysisOpt <- foreach(kk = 1:length(AllDrugPredictModelsNew)) %dopar% {
  # scaled data
  rowPosition<-which(AllDrugPredictModelsNew[[kk]]$scale$opt.param$alpha.opt == alpha)
  colPosition<-which(AllDrugPredictModelsNew[[kk]]$scale$opt.param$lambda.opt == lambda)
  
  mseTrain<-c()
  mseTest<-c()
  for (k in 1:numFolds){
    mseTrain<-c(mseTrain,mean((AllDrugPredictModelsNew[[kk]]$scale$cv.all[[rowPosition]]$foldTrainPredictions[[k]][,colPosition] 
                               - AllDrugPredictModelsNew[[kk]]$scale$cv.all[[rowPosition]]$foldTrainObservations[[k]])^2))
    mseTest<-c(mseTest,mean((AllDrugPredictModelsNew[[kk]]$scale$cv.all[[rowPosition]]$foldTestPredictions[[k]][,colPosition] 
                             - AllDrugPredictModelsNew[[kk]]$scale$cv.all[[rowPosition]]$foldTestObservations[[k]])^2))
  }
  # unscaled data
  rowPosition<-which(AllDrugPredictModelsNew[[kk]]$nonscale$opt.param$alpha.opt == alpha)
  colPosition<-which(AllDrugPredictModelsNew[[kk]]$nonscale$opt.param$lambda.opt == lambda)
  
  mseTrain1<-c()
  mseTest1<-c()
  for (k in 1:numFolds){
    mseTrain1<-c(mseTrain1,mean((AllDrugPredictModelsNew[[kk]]$nonscale$cv.all[[rowPosition]]$foldTrainPredictions[[k]][,colPosition] 
                               - AllDrugPredictModelsNew[[kk]]$nonscale$cv.all[[rowPosition]]$foldTrainObservations[[k]])^2))
    mseTest1<-c(mseTest1,mean((AllDrugPredictModelsNew[[kk]]$nonscale$cv.all[[rowPosition]]$foldTestPredictions[[k]][,colPosition] 
                             - AllDrugPredictModelsNew[[kk]]$nonscale$cv.all[[rowPosition]]$foldTestObservations[[k]])^2))
  }
  return(list(scale = cbind(mseTrain,mseTest),nonscale = cbind(mseTrain1,mseTest1)))
}
