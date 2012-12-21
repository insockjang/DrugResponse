require(SGL)
mySGLModel <- setRefClass(Class = "mySGLModel",
                          contains="PredictiveModel",
                          fields="model",
                          methods = list(
                            initialize = function(...){
                              return(.self)
                            },
                            
                            rawModel = function(){
                              return(.self$model)
                            },
                            
                            customTrain = function(featureData, responseData, index = index,nfolds,...){
                              
                              .self$model<-cvSGL(list(y=responseData,x=featureData),index = index,type = "linear",alpha = 0,nfold = nfolds, ...)
                              .self$model$optLambda<-(.self$model$lambdas[which.min(.self$model$lldiff)])
                              .self$model$coeff<-.self$model$fit$beta[,which.min(.self$model$lldiff)]
                              
                            },
                            
                            customPredict = function(featureData){
                              predictedResponse <-  (featureData %*% .self$model$coeff)
                              return(predictedResponse)
                            },
                            
                            getCoefficients = function(){
                              .self$model$coeff<-.self$model$fit$beta[,which.min(.self$model$lldiff)]
                              return(.self$model$coeff)
                            }
                            
                            )
                          
                          )
