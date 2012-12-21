#' Perform cross validation on a preditive model, this model support Cox model and predictiveModeling
#'
#' @param featureData
#' @param responseData
#' @param model either an instance of PredictiveModel or a string holding the name of one of the machine learning methods that caret supports
#' @param numFolds defaults to 5
#' @param trControl defaults to defaultTrainControl
#' @return a list of PredictiveModelPerformance one per fold
#' @seealso defaultTrainControl
#' @export
crossValidateFindOptimalCoxModel <- 
  function(featureData, responseData, model, numFolds = 5, GRID = NULL, ...){
    
    if (is(model, "character")){
      model <- CaretModel$new(modelType = model)
    }
    
    #-----------------------------------------------------------------------
    # Split the data into training and test partitions
    # -----------------------------------------------------------------------
    # set.seed(2)
    foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
    
    message("Training partitions cross validation step in optimization")
    
    modelResults <- foreach(fold = foldIndices) %dopar% {
      foldModel <- model$copy()
      foldModel$customTrain(featureData[-fold,], responseData[-fold], alpha = GRID[1],lambda = GRID[2])
      return(foldModel)
    }
    
    foldTestPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[foldIndices[[i]],])
    #testPredictions <- do.call("rbind", foldTestPredictions)
    
    foldTrainPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[-foldIndices[[i]],])
    #trainPredictions <- do.call("rbind", foldTrainPredictions)
    
    foldTestObservations <- foreach(i=1:numFolds) %do% responseData[foldIndices[[i]]]
    #testObservations <- do.call("rbind", foldTestObservations)
    
    foldTrainObservations <- foreach(i=1:numFolds) %do% responseData[-foldIndices[[i]]]
    #trainObservations <- do.call("rbind", foldTrainObservations)
    
    res <- list(trainPredictions = foldTrainPredictions, trainObservations = foldTrainObservations,
                testPredictions = foldTestPredictions, testObservations = foldTestObservations)
    
    return(res)
  }
