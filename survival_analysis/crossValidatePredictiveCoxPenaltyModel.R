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
crossValidatePredictiveCoxPenaltyModel <- 
  function(featureData, responseData, model, numFolds = 5, trControl = defaultTrainControl(), penaltyFactor = rep(1,nvars),...){
    
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
      foldModel$customTrain(featureData[-fold,], responseData[-fold], penaltyFactor = penaltyFactor)
      return(foldModel)
    }
    coefficients <- foreach(i=1:numFolds) %do% coef(modelResults[[i]]$rawCaretModel(),s="lambda.min")
    
    foldTestPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[foldIndices[[i]],])
    #testPredictions <- do.call("rbind", foldTestPredictions)
    
    foldTrainPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[-foldIndices[[i]],])
    #trainPredictions <- do.call("rbind", foldTrainPredictions)
    
    foldTestObservations <- foreach(i=1:numFolds) %do% responseData[foldIndices[[i]]]
    #testObservations <- do.call("rbind", foldTestObservations)
    
    foldTrainObservations <- foreach(i=1:numFolds) %do% responseData[-foldIndices[[i]]]
    #trainObservations <- do.call("rbind", foldTrainObservations)
    
    res <- list(trainPredictions = foldTrainPredictions, trainObservations = foldTrainObservations,
                testPredictions = foldTestPredictions, testObservations = foldTestObservations, coefficients = coefficients)
    
    return(res)
  }
