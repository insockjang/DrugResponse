# for testing customized Class for partial least square as one testbed for predictiveModeling
setClass(
  Class = 'MyPLSModel',
  contains="PredictiveModel",
  representation = representation(
    model='list')
)

# method is already implemented and utilized(wrapped up) in caret
# This is for training
setMethod(
  f = "customTrain",
  
  signature = signature("MyPLSModel", "matrix", "numeric"),
 
  definition = function(method, featureData, responseData){
    
    message("in customTrain ", class(method))

    nCompGrid <- expand.grid(.ncomp = c(1,2,3,4))
    caretTrainModel <- caret::train(featureData, 
                          responseData, 
                          method = "pls",
                          preProcess = NULL, # preProcess = c("center", "scale"),
                          #tuneLength = 4,
                          trControl = defaultTrainControl(),
                          tuneGrid = nCompGrid
                          )
    
    modelHolder <- list()
    modelHolder$model <- caretTrainModel
    
    method@model <- modelHolder
    
    method
  }
)

# This is for prediction(testing)
setMethod(
  
  f = "customPredict",
  
  signature = signature("MyPLSModel", "matrix"),
  
  definition = function(method, featureData){
  
    message("in customPredict ", class(method))

    predict(method@model$model, featureData)
  }
)
