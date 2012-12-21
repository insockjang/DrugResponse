cvGridAlphaMyENetModel<-function(featureData, responseData, model, alphas = alphas, lambda = lambda, numFolds = 3, ...){
    
    foldIndices = createFolds(featureData[,1],k = numFolds)
  
    foldResults <-foreach(fold = foldIndices) %dopar% {
      return(gridAlphaENetResults(featureData[-fold,],responseData[-fold],model = myEnetModel,alphas = alphas,lambda = lambda, nfolds =5))
    }  
    
    foldTestPredictions <- foreach(i = 1:numFolds) %do% foldResults[[i]]$bestModel$customPredict(featureData[foldIndices[[i]],])
    testPredictions <- do.call("c", foldTestPredictions)
    
    foldTestObservations <- foreach(i = 1:numFolds) %do% responseData[foldIndices[[i]]]
    testObservations <- do.call("c", foldTestObservations)
    
    foldTrainPredictions <- foreach(i = 1:numFolds) %do% foldResults[[i]]$bestModel$customPredict(featureData[-foldIndices[[i]],])
    trainPredictions <- do.call("c", foldTrainPredictions)
    
    foldTrainObservations <- foreach(i = 1:numFolds) %do% responseData[-foldIndices[[i]]]
    trainObservations <- do.call("c", foldTrainObservations)
    
    res <- list(trainPredictions = trainPredictions, trainObservations = trainObservations,
                testPredictions = testPredictions, testObservations = testObservations)
    
    return(res)
}
