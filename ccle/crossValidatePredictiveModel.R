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
crossValidatePredictiveModel <- 
  function(featureData, responseData, model, numFolds = 5, trControl = defaultTrainControl(), ...){
    
    if (is(model, "character")){
      model <- CaretModel$new(modelType = model)
    }
    
    #-----------------------------------------------------------------------
    # Split the data into training and test partitions
    # -----------------------------------------------------------------------
    set.seed(2)
    foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
    
    message("Training partitions")
    
    modelResults <- foreach(fold = foldIndices) %dopar% {
      foldModel <- model$copy()
      foldModel$customTrain(featureData[-fold,], responseData[-fold])
      return(foldModel)
    }
    
    foldTestPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[foldIndices[[i]],])
    
    foldTrainPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[-foldIndices[[i]],])
    
    foldTestObservations <- foreach(i=1:numFolds) %do% responseData[foldIndices[[i]]]
    
    foldTrainObservations <- foreach(i=1:numFolds) %do% responseData[-foldIndices[[i]]]
    
    res <- list(trainPredictions = foldTrainPredictions, trainObservations = foldTrainObservations,
                testPredictions = foldTestPredictions, testObservations = foldTestObservations)
    
    return(res)
  }
