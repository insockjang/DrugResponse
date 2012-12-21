#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
require(glmnet)
require(survival)
myEnetCoxModel <- setRefClass(Class = "myEnetCoxModel",
                          contains="PredictiveModel",
                          fields="model",
                          methods = list(
                            initialize = function(...){
                              return(.self)
                            },
                            
                            rawModel = function(){
                              return(.self$model)
                            },
                            
                            customTrain = function(featureData, responseData, alpha, lambda, nfolds,...){
                              lambda = sort(lambda,decreasing = TRUE)
                              foldIndices = createFolds(featureData[,1],k = nfolds)
                              if(length(alpha) == 1){
                                results<-c()
                                for(fold in foldIndices){
                                  fit<-glmnet(featureData[-fold,],responseData[-fold,],family = "cox", alpha = alpha, lambda = lambda,...) 
                                  pred<-predict(fit,featureData[fold,])
                                  cIndex<-c()
                                  for(k in 1:ncol(pred)){
                                    cIndex<-c(cIndex,survConcordance(responseData[fold,]~pred[,k])$concordance)
                                  }
                                  results <- rbind(results, cIndex)
                                }
                                metrics <- apply(results,2,mean)
                                optParam <- c(max(metrics), alpha, lambda[which.max(metrics)])
                                names(optParam) <- c("cIndex","alpha","lambda")
                                .self$model <- glmnet(featureData,responseData,family = "cox", alpha = optParam[2], lambda = optParam[3],...)
                                .self$model$optParam <- optParam                                          
                                
                              }
                              else{
                                optParam <-c()
                                for(kk in 1:length(alpha)){                                  
                                  results<-c()
                                  for(fold in foldIndices){
                                    fit<-glmnet(featureData[-fold,],responseData[-fold,],family = "cox", alpha = alpha[kk], lambda = lambda,...) 
                                    pred<-predict(fit,featureData[fold,])
                                    cIndex<-c()
                                    for(k in 1:ncol(pred)){
                                      cIndex<-c(cIndex,survConcordance(responseData[fold,]~pred[,k])$concordance)
                                    }
                                    results<- rbind(results, cIndex)
                                  }
                                  metrics <- apply(results,2,mean)
                                  vec<-c(max(metrics),alpha[kk],lambda[which.max(metrics)])
                                  optParam<-rbind(optParam,vec)
                                }
                                colnames(optParam) <- c("cIndex","alpha","lambda")
                                bestModel <-which.max(optParam[,1])
                                .self$model <- glmnet(featureData,responseData,family = "cox", alpha = optParam[bestModel,2],lambda = optParam[bestModel,3],...)
                                .self$model$optParam = optParam
                                
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
