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
bootstrapPredictiveCoxModel <- 
  function(featureData, responseData, model, numBoots = 100, trControl = defaultTrainControl(),...){
    
    if (is(model, "character")){
      model <- CaretModel$new(modelType = model)
    }
    
    #-----------------------------------------------------------------------
    # Resample the data
    # -----------------------------------------------------------------------
    set.seed(2)
    bootIndices <- createResample(featureData[,1], times = numBoots, list = TRUE)
    
    Coefs <- foreach(boot = bootIndices) %dopar% {
      bootModel <- model$copy()
      bootModel$customTrain(featureData[boot,], responseData[boot])
      return(as.logical(coef(bootModel$rawCaretModel(),s="lambda.min")))      
    }
    CoefResults <- do.call("cbind", Coefs)
    return(CoefResults)
  }
