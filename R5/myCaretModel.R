#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
myCaretModel <- setRefClass(Class = "myCaretModel",
                          contains="PredictiveModel",
                          fields=c("model", "modelType"),
                          methods = list(
                            initialize = function(modelType, ...){
                              .self$modelType = modelType
                              
                              return(.self)
                            },
                            
                            rawCaretModel = function(){
                              return(.self$model)
                            },
                            
                            copy = function() {
                              
                              result <- CaretModel$new(.self$modelType)
                              result$model <- .self$model
                              
                              return(result)
                            },
                            
                            customTrain = function(featureData, responseData, trControl, tuneGrid, ...){
                              .self$model <- caret::train(featureData, 
                                                          responseData, 
                                                          method = .self$modelType,
                                                          preProcess = NULL, # preProcess = c("center", "scale"),
                                                          tuneLength = 4,
                                                          trControl = trControl,
                                                          tuneGrid = tuneGrid,
                                                          ...
                                                          )
                            },
                            
                            customPredict = function(featureData){
                              predictedResponse <- predict(.self$model, featureData)
                              return(predictedResponse)
                            }
                            )
                          )
