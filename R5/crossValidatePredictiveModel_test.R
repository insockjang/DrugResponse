crossValidatePredictiveModel_test <- 
  function(featureData, responseData, model, numFolds = 5, ...){
    
    #-----------------------------------------------------------------------
    # Split the data into training and test partitions
    # -----------------------------------------------------------------------
    set.seed(2)
    foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
    
    message("Training partitions")
    
    foldResults <- foreach(fold = foldIndices) %dopar% {   
      
      model$customTrain(featureData[-fold,], responseData[-fold], ...)
      
      trainPredictions = model$customPredict(featureData[-fold,])
      trainObservations = responseData[-fold]
      testPredictions = model$customPredict(featureData[fold,])
      testObservations = responseData[fold]
      
      res <- list(trainPredictions = trainPredictions, 
                  trainObservations = trainObservations,
                  testPredictions = testPredictions,
                  testObservations = testObservations)
      return(res)           
    }    
    return(foldResults)
  }
