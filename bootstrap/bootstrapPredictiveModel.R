#' Perform bootstrapping on a preditive model in order to select more powerful features, this method support both caret models and PredictiveModels
#'
#' @param featureData
#' @param responseData
#' @param method either an instance of PredictiveModel or a string holding the name of one of the machine learning methods that caret supports
#' @param numBootstrap defaults to 100
#' @return a vector of features which shows how many times selected 
#' @export
bootstrapPredictiveModel <- 
  function(featureData, responseData, model, numBootstrap = 100, filterData = TRUE, trControl = trControl, ...){
    if (filterData == TRUE){
      processedData <- filterPredictiveModelData(featureData, responseData)
      featureData <- processedData$featureData
      responseData <- processedData$responseData
    }
    
    bootIndices <- createResample(featureData[,1],times = numBootstrap,list = TRUE)
    bootResults <- foreach (numBoot = bootIndices) %dopar% {
      model$customTrain(featureData[numBoot, ], responseData[numBoot], trControl = trControl,...)
      caretModel <- model$rawCaretModel()
      
      # coef_vecBootstrap is for binding column 
      # coef_binBootstrap is for counting cumulate selected features in each bootstrapping step
      return(list(coef_vecBootstrap = caretModel$finalModel$beta[,ncol(caretModel$finalModel$beta)],
                  coef_binBootstrap = caretModel$finalModel$beta[,ncol(caretModel$finalModel$beta)]!=0))
    }
    
    bootCoef <- foreach(i = 1:numBootstrap) %do% bootResults[[i]]$coef_binBootstrap
    cumBootCoef <- do.call("cbind", bootCoef)
    sumBootCoef <-apply(cumBootCoef,1,sum)
    return(sumBootCoef)
  }
