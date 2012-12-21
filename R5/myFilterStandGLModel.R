require(standGL)
myFilterStandGLModel <- setRefClass(Class = "myFilterStandGLModel",
                                    contains="PredictiveModel",
                                    fields=c("beta","param"),
                                    methods = list(
                                      initialize = function(...){
                                        return(.self)
                                      },
                                      
                                      rawModel = function(){
                                        return(.self$model)
                                      },
                                      
                                      customTrain = function(featureData, responseData, index = index, alpha, lambda, nfolds,...){
                                        lambda = sort(lambda,decreasing = TRUE)
                                        foldIndices = createFolds(featureData[,1],k = nfolds,list = TRUE)
                                        
                                        MSE<-foreach(fold = foldIndices) %dopar% {                                      
                                          fData<-filterPredictiveModelData(featureData[-fold,],responseData[-fold])
                                          FeatureData<-fData$featureData
                                          ResponseData<-fData$responseData
                                          
                                          fit<-standGL(ResponseData,FeatureData,index = index,alpha = 1,lam.path = lambda)
                                          
                                          pred <-t(fit$beta) %*% t(featureData[fold,colnames(FeatureData)])
                                          
                                          rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
                                          
                                          mse<-c()
                                          for(k in 1:dim(pred)[1]){
                                            mse<-c(mse,rmse(pred[k,],responseData[fold]))
                                          }
                                          return(mse)
                                        }
                                        mse <- do.call("rbind", MSE)
                                        metrics<-apply(mse,2,mean)                                    
                                        optParam <- c(min(metrics), alpha, lambda[which.min(metrics)])
                                        
                                        names(optParam) <- c("RMSE","alpha","lambda")
                                        FIT <-standGL(responseData,featureData,index = index,alpha = 1,lam.path = lambda)
                                        .self$beta <- FIT$beta[,which(lambda == optParam[3])]
                                        .self$param <- optParam
                                        
                                      },
                                      
                                      customPredict = function(featureData){
                                        predictedResponse <-  (featureData %*% .self$beta)
                                        return(predictedResponse)
                                      },
                                      
                                      getCoefficients = function(){
                                        return(.self$beta)
                                      }
                                      
                                      )
                                    
                                    )
