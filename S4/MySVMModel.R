setClass(
  Class = 'MySVMModel',
  contains="PredictiveModel",
  representation = representation(
    model='list')
  )

setMethod(
  f = "customTrain",
  signature = signature("MySVMModel", "matrix", "numeric"),
  definition = function(method, featureData, responseData){
    message("in customTrain ", class(method))
    btControl <- trainControl(number = floor(ncol(featureData)*0.2))
    caretTrainModel <- caret::train(featureData, 
                                    responseData, 
                                    method = "svmRadial",
                                    tuneLength = 5,
                                    trControl = btControl, # "number" means No. of bootstrapping
                                    scaled = FALSE
                                    )
    
    modelHolder <- list()
    modelHolder$model <- caretTrainModel
    
    method@model <- modelHolder
    
    method
  }
  )

setMethod(
  f = "customPredict",
  signature = signature("MySVMModel", "matrix"),
  definition = function(method, featureData){
    message("in customPredict ", class(method))
    predict(method@model$model, featureData)
  }
  )
