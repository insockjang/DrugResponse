#' Perform bootstrapping on a preditive model in order to select more powerful features, this method support both caret models and PredictiveModels
#'
#' @param featureData
#' @param responseData
#' @param method either an instance of PredictiveModel or a string holding the name of one of the machine learning methods that caret supports
#' @param numBootstrap defaults to 100
#' @return a vector of features which shows how many times selected 
#' @export
bootstrapPredictiveModel <- 
  function(featureData, responseData, model, numBootstrap = 100, ...){
    
    bootIndices <- createResample(featureData[,1],times = numBootstrap,list = TRUE)
    
    bootResults <- foreach (numBoot = bootIndices) %dopar% {
      bootModel<-model$copy()
      bootModel$customTrain(featureData[numBoot, ], responseData[numBoot], ...)
      coeffs<-bootModel$getCoefficients()      
      return(coeffs)
    }
    
    return(bootResults)
  }
