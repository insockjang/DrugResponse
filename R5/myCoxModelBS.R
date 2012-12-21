#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
require(glmnet)

myCoxModelBS <- setRefClass(Class = "myCoxModelBS",
                          contains="PredictiveModel",
                          fields="coefficientBS",
                          methods = list(
                            initialize = function(...){
                              return(.self)
                            },
                            
                            rawModel = function(){
                              return(.self$model)
                            },
                            
                            customTrain = function(featureData, responseData, alpha, lambda, numBootstrap,...){
                              lambda = sort(lambda,decreasing = TRUE)
                              
                              bootIndices <- createResample(featureData[,1],times = numBootstrap,list = TRUE)
                              
                              if(length(alpha) == 1){
                                coefficients<-list()
                                for(kk in 1:numBootstrap){
                                  fit<-cv.glmnet(featureData[bootIndices[[kk]],],responseData[bootIndices[[kk]],],family = "cox",type.measure="deviance",alpha = alpha,lambda = lambda, grouped= TRUE)                                   
                                  coefficients[[kk]]<-coef(fit,s= "lambda.min")
                                }
                                .self$coefficientBS <- coefficients
                                                                                                      
                              }
                              else{
                               coefficients <-list()
                               for(kk in 1:numBootstrap){
                                 optParam <-c()
                                 for(kk in 1:length(alpha)){                                                                  
                                   fit<-cv.glmnet(featureData[bootIndices[[kk]],],responseData[bootIndices[[kk]],],family = "cox", type.measure="deviance",alpha = alpha[kk],lambda = lambda, grouped= TRUE)                                   
                                   metrics <- fit$cvm
                                   vec<-c(min(metrics),alpha[kk],lambda[which.min(metrics)])
                                   optParam<-rbind(optParam,vec)
                                 }
                                 colnames(optParam) <- c("Partial Likelihood Deviance","alpha","lambda")
                                 bestModel <-which.min(optParam[,1])
                                 
                                 fit <- glmnet(featureData,responseData,family = "cox", alpha = optParam[bestModel,2],lambda = optParam[bestModel,3],...)
                                 coefficients[[kk]]<-coef(fit)
                               }
                                .self$coefficientBS = coefficients
                                
                              }
                            },
                            
                            customPredict = function(){
                              K=length(.self$coefficientBS)
                              
                              featureName<-rownames(.self$coefficientBS[[1]])
                              A<-matrix(0,nrow = length(featureName), ncol = K)
                              for(k in 1:K){
                                A[,k] = as.numeric(.self$coefficientBS[[k]])
#                                 A[,k] = as.numeric(fit$coefficientBS[[k]])
                              }
                              B<- A!=0
                              BS<-apply(B,1,sum)
                              names(BS)<-featureName
                              return(BS)                              
                            }
                            
                            )
                          )
