filterPredictiveModelData1<-function(dataSet,kk,refineFilter = TRUE){
  
  if(refineFilter){
    filteredData <- filterPredictiveModelData(dataSet$featureData,dataSet$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  }
  else{
    filteredData <- filterPredictiveModelData(dataSet$featureData,dataSet$responseData[,kk,drop=FALSE])
  }
  
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData  
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)
  return(list(featureData = filteredFeatureDataScaled, responseData = filteredResponseDataScaled))
}
