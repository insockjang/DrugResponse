### DEMO ICR
library(predictiveModeling)
library(synapseClient)
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


AllDrugPredictModelIcr<-foreach(kk = 1:ncol(responseDataValidation)) %dopar% {
  
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
  
  
  # 10 fold cross validation to find optimal number of principal components    
  myTrControl=trainControl(method = "cv", number = 5, returnResamp = "all", verboseIter = TRUE)
  
  icrScaled <- CaretModel$new(modelType = "icr")
  icrUnscaled <- CaretModel$new(modelType = "icr")
  
  icrScaled$customTrain(filteredFeatureDataTrain_scaled, filteredResponseDataTrain_scaled, tuneGrid=expand.grid(.n.comp = seq(1,100,1)), trControl = myTrControl)
  predTrain <- icrScaled$customPredict(filteredFeatureDataTrain_scaled)
  obsTrain  <- filteredResponseDataTrain_scaled
  predTest <- icrScaled$customPredict(filteredFeatureDataValidation_scaled)
  obsTest  <- filteredResponseDataValidation_scaled
  
  icrUnscaled$customTrain(filteredFeatureDataTrain, filteredResponseDataTrain, tuneGrid=expand.grid(.n.comp = seq(1,100,1)), trControl = myTrControl)
  pred1Train <- icrUnscaled$customPredict(filteredFeatureDataTrain)
  obs1Train  <- filteredResponseDataTrain
  pred1Test <- icrUnscaled$customPredict(filteredFeatureDataValidation)
  obs1Test  <- filteredResponseDataValidation
  
  return(list(ScaledModel = icrScaled, 
              UnscaledModel = icrUnscaled, 
              testPredictionsScaled = predTest,
              testObservationsScaled = obsTest,
              testPredictionsUnscaled = pred1Test,
              testObservationsUnscaled = obs1Test,
              trainPredictionsScaled = predTrain,
              trainObservationsScaled = obsTrain,
              trainPredictionsUnscaled = pred1Train,
              trainObservationsUnscaled = obs1Train
              )
         )    
}



#########################################################################################################
########  Analysis Step: Training and Testing data are scaled(normalized) vs. raw(unnormalized)  ########
#########################################################################################################
# AllDrugPredictModelPcr 1~24 should be checked with following code.
correlationScaled<-c()
correlationUnscaled<-c()

for(k in 1:length(AllDrugPredictModelIcr)){
  correlationScaled<-c(correlationScaled,cor(AllDrugPredictModelIcr[[k]]$foldScaledPredictions,AllDrugPredictModelIcr[[k]]$foldScaledObservations, use = "na.or.complete"))
  correlationUnscaled<-c(correlationUnscaled,cor(AllDrugPredictModelIcr[[k]]$foldUnscaledPredictions,AllDrugPredictModelIcr[[k]]$foldUnscaledObservations, use = "na.or.complete"))
}
