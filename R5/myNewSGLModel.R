require(SGL)
myNewSGLModel <- setRefClass(Class = "myNewSGLModel",
                                    contains="PredictiveModel",
                                    fields="model",
                                    methods = list(
                                      initialize = function(...){
                                        return(.self)
                                      },
                                      
                                      rawModel = function(){
                                        return(.self$model)
                                      },
                                      
                                      customTrain = function(featureData, responseData, index = index, alpha, nfolds,...){                                      
                                        
                                        sgl<-foreach(k = 1:length(alphas)) %dopar% {
                                          fit<-cvSGL(list(y=responseData,x=featureData),index = index,type = "linear",alpha = alphas[k],nfold = nfolds)
                                          return(fit)
                                        }
                                        
                                        llDiff<-foreach(k = 1:length(alpha)) %do% {sgl[[k]]$lldiff}
                                        LL<-do.call("rbind",llDiff)
                                        
                                        .self$model<-sgl[[which.min(apply(LL,1,min))]]
                                        .self$model$optAlpha<-alpha[which.min(apply(LL,1,min))]
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
