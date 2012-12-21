#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @export
StandGLModel <- setRefClass(Class = "StandGLModel",
                              contains="PredictiveModel",
                              fields="model",
                              methods = list(
                                initialize = function(...){
                                  return(.self)
                                },
                                
                                rawModel = function(){
                                  return(.self$model)
                                },
                                
                                customTrain = function(featureData, responseData, index = index, alpha = alpha, lambda = lambda, ...){
                                  .self$model <- standGL(responseData, featureData, index = index, alpha = alpha, lam.path = lambda)                                   
                                },
                                
                                customPredict = function(featureData){
                                  predictedResponse <-  featureData %*% .self$model$beta
                                  return(predictedResponse)
                                },
                                
                                getCoefficients = function(){
                                  return(.self$model$beta)
                                }
                                
                                )
                              
                              )
