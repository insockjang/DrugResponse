#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @export
myCVStandGLModel <- setRefClass(Class = "myCVStandGLModel",
                              contains="PredictiveModel",
                              fields="model",
                              methods = list(
                                initialize = function(...){
                                  return(.self)
                                },
                                
                                rawModel = function(){
                                  return(.self$model)
                                },
                                
                                customTrain = function(featureData, responseData, index = index, lambda = lambda, ...){
                                #customTrain = function(featureData, responseData, index = index, lambda = lambda, ...){
                                    #intercept<-rep(1,length(responseData))
                                  #iFeatureData<-cbind(intercept,featureData)
                                  #iIndex <- c(0,index) +1
                                  #iIs.pen <- c(0,rep(1,ncol(featureData)))
                                  #fit <- cv.standGL(responseData, iFeatureData, index = iIndex, alpha = 1, lam.path = lambda, is.pen = iIs.pen, nfold = 3) 
                                  .self$model <- cv.standGL(responseData, featureData, index = index, alpha = 1, lam.path = lambda, nfold = 3) 
                                  #.self$model <-list(intercept=fit$beta[1,],beta = fit$beta[-1,],lambda = fit$lam.path)
                                },
                                
                                customPredict = function(featureData){
                                  predictedResponse <-  cbind(rep(1,ncol(featureData)),featureData) %*% rbind(.self$model$intercept,.self$model$beta)
                                  return(predictedResponse)
                                },
                                
                                getCoefficients = function(){
                                  return(.self$model$beta)
                                }
                                
                                )
                              
                              )
