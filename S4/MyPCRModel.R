setClass(
  Class = 'MyPCRModel',
  contains="PredictiveModel",
  representation = representation(
    model='list')
)

setMethod(
  
  f = "customTrain",
  
  signature = signature("MyPCRModel", "matrix", "numeric"),
  
  definition = function(method, featureData, responseData){

    message("in customTrain ", class(method))
    
    nCompGrid <- expand.grid(.ncomp = c(1,2,3,4))
    
    caretTrainModel <- caret::train(featureData, 
                          responseData, 
                          method = "pcr",
                          preProcess = NULL, # preProcess = c("center", "scale"),
                          #tuneLength = 4,
                          tuneGrid = nCompGrid
                          )
    
    modelHolder <- list()
    modelHolder$model <- caretTrainModel
    
    method@model <- modelHolder
    
    method
  }
  )

setMethod(
  
  f = "customPredict",
  
  signature = signature("MyPCRModel", "matrix"),
  
  definition = function(method, featureData){
    
    message("in customPredict ", class(method))

    predict(method@model$model, featureData)
  }
)
