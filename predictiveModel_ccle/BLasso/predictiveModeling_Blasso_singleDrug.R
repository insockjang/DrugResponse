### DEMO Blasso
# Be careful : Bayesia Lasso/Ridge need not crossvalidation. 
# CVs are needed to decide optimal parameters.
library(predictiveModeling)
library(synapseClient)
library(monomvn)
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myBlassoModel.R")
synapseLogin("in.sock.jang@sagebase.org")

id_Layer <- "158966"     
layer <- loadEntity(id_Layer)
layer$objects$trainData



# 5 times of Model validation with randomly splitted dataset
featureDataTrain <- createAggregateFeatureDataSet(list(expr = layer$objects$trainData$exprData, 
                                                       copy = layer$objects$trainData$copyData, 
                                                       mut = layer$objects$trainData$oncomapData))
responseDataTrain <- layer$objects$trainData$drugData


featureDataValidation <- createAggregateFeatureDataSet(list(expr = layer$objects$testData$exprData, 
                                                            copy = layer$objects$testData$copyData, 
                                                            mut = layer$objects$testData$oncomapData))
responseDataValidation <- layer$objects$testData$drugData


dataSets_ccle <- createFeatureAndResponseDataList(t(featureDataTrain),responseDataTrain)

for(kk in 1:ncol(responseDataValidation)){  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredDataTrain<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # training data
  filteredFeatureDataTrain  <- filteredDataTrain$featureData
  filteredResponseDataTrain <- filteredDataTrain$responseData
  
  # validating data
  filteredFeatureDataValidation <- t(featureDataValidation[colnames(filteredFeatureDataTrain),])
  filteredResponseDataValidation <- responseDataValidation[,kk][[1]]
  
  ## scale this data    
  filteredFeatureDataTrain_scaled <- scale(filteredFeatureDataTrain)
  filteredResponseDataTrain_scaled <- scale(filteredResponseDataTrain)  
  
  ## scale validation data    
  filteredFeatureDataValidation_scaled <- scale(filteredFeatureDataValidation)
  filteredResponseDataValidation_scaled <- scale(filteredResponseDataValidation)  
  
  
  blassoScaled <- myBlassoModel$new()
  blassoUnscaled <- myBlassoModel$new()
  
  blassoScaled$customTrain(filteredFeatureDataTrain_scaled, filteredResponseDataTrain_scaled, mcmcNum = 10000, burnIn = 2000)
  predTrain <- blassoScaled$customPredict(filteredFeatureDataTrain_scaled)
  predTest <- blassoScaled$customPredict(filteredFeatureDataValidation_scaled)
  obsTrain  <- filteredResponseDataTrain_scaled
  obsTest  <- filteredResponseDataValidation_scaled
  
  filename <- paste("drug_",kk,"_BlassoScaled_10000.Rdata",sep ="")
  
  results <-c()
  
  results<-list(ScaledModel = blassoScaled, 
                testPredictionsScaled = predTest,
                testObservationsScaled = obsTest,
                trainPredictionsScaled = predTrain,
                trainObservationsScaled = obsTrain
                )
  save(results,file = filename)
  
  blassoUnscaled$customTrain(filteredFeatureDataTrain, filteredResponseDataTrain, mcmcNum = 10000, burnIn = 2000)
  pred1Train <- blassoUnscaled$customPredict(filteredFeatureDataTrain)
  pred1Test <- blassoUnscaled$customPredict(filteredFeatureDataValidation)
  obs1Train  <- filteredResponseDataTrain
  obs1Test  <- filteredResponseDataValidation
  
  filename <- paste("drug_",kk,"_BlassoUnscaled_10000.Rdata",sep ="")
  
  results <-c()
  
  results<-list(UnscaledModel = blassoUnscaled, 
                testPredictionsUnscaled = pred1Test,
                testObservationsUnscaled = obs1Test,
                trainPredictionsUnscaled = pred1Train,
                trainObservationsUnscaled = obs1Train
                )
  save(results,file = filename)
}




#########################################################################################################
########  Analysis Step: Training and Testing data are scaled(normalized) vs. raw(unnormalized)  ########
#########################################################################################################
# AllDrugPredictModelPcr 1~24 should be checked with following code.
correlationScaled<-c()
correlationUnscaled<-c()

for(k in 1:length(AllDrugPredictModelIcr)){
  correlationScaled<-c(correlationScaled,cor(AllDrugPredictModelblasso[[k]]$foldScaledPredictions,AllDrugPredictModelblasso[[k]]$foldScaledObservations, use = "na.or.complete"))
  correlationUnscaled<-c(correlationUnscaled,cor(AllDrugPredictModelblasso[[k]]$foldUnscaledPredictions,AllDrugPredictModelblasso[[k]]$foldUnscaledObservations, use = "na.or.complete"))
}
