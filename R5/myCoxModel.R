#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
require(glmnet)

myCoxModel <- setRefClass(Class = "myCoxModel",
                        contains="PredictiveModel",
                        fields="model",
                        methods = list(
                             initialize = function(...){
                               return(.self)
                             },
                             
                             rawModel = function(){
                               return(.self$model)
                             },
                             
                             customTrain = function(featureData, responseData, alpha, lambda,...){
                               lambda = sort(lambda,decreasing = TRUE)
                               if(length(alpha) == 1){
                                 fit<-cv.glmnet(featureData,responseData,family = "cox",type.measure="deviance",alpha = alpha,lambda = lambda, grouped= TRUE)                                   
                                 metrics <- fit$cvm
                                 optParam <- c(min(metrics), alpha, lambda[which.min(metrics)])
                                 names(optParam) <- c("Partial Likelihood Deviance","alpha","lambda")
                                 .self$model <- glmnet(featureData,responseData,family = "cox", alpha = optParam[2], lambda = optParam[3],...)
                                 .self$model$optParam <- optParam                                                                         
                               }
                               else{
                                 optParam <-c()
                                 for(kk in 1:length(alpha)){                                                                  
                                   fit<-cv.glmnet(featureData,responseData,family = "cox", type.measure="deviance",alpha = alpha[kk],lambda = lambda, grouped= TRUE)                                   
                                   metrics <- fit$cvm
                                   vec<-c(min(metrics),alpha[kk],lambda[which.min(metrics)])
                                   optParam<-rbind(optParam,vec)
                                 }
                                 colnames(optParam) <- c("Partial Likelihood Deviance","alpha","lambda")
                                 bestModel <-which.min(optParam[,1])
                                 .self$model <- glmnet(featureData,responseData,family = "cox", alpha = optParam[bestModel,2],lambda = optParam[bestModel,3],...)
                                 .self$model$Param = optParam
                                 .self$model$optParm = optParam[bestModel,]
                                 
                               }
                             },
                             customPredict = function(featureData){
                               predictedResponse <- predict(.self$model,featureData)
                               return(predictedResponse)
                             },
                             
                             getCoefficients = function(){
                               return(coef(.self$model))                              
                             }
                             
                             
                             )
                           )
