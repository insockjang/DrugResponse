#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
ENetModel <- setRefClass(Class = "ENetModel",
                           contains="PredictiveModel",
                           fields="model",
                           methods = list(
                             initialize = function(...){
                               return(.self)
                             },
                             
                             rawModel = function(){
                               return(.self$model)
                             },
                             customTrain = function(featureData, responseData, alpha, lambda, ...){
                               .self$model <- glmnet(featureData,responseData, alpha = alpha, lambda = lambda, ...)
                             },
                             
                             customPredict = function(featureData){
                               predictedResponse <- predict(.self$model, featureData)
                               return(predictedResponse)
                             },
                             
                             getCoefficients = function(){
                               return(coef(.self$model))
                             }
                             
                             )
                           
                           )
