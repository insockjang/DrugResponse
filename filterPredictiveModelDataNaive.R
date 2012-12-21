#' filterPredictiveModelDataNaive filters ... NEED MORE HERE
#'
#' @param featureData
#' @param responseData
#' @return a list of the filtered feature and response data
#' @export
filterPredictiveModelDataNaive <- 
  function(featureData, responseData)
  {
    isNas <- is.na(responseData)
    
    featureData <- featureData[!isNas,]
    responseData <- responseData[!isNas]
    
    variances <- apply(featureData, 2, var)
    if(any(variances==0 | is.na(variances) )){
      featureData <- featureData[, -which(variances==0 | is.na(variances))] 
    }
    list(featureData = featureData, responseData = responseData)
  }
