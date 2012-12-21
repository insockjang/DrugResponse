#' Constructor for a class in the PredictiveModel class hierarchy that implements a simple function regressing the top correlated features
#' R5 Class for Bayesian Ridge for predictiveModeling
#' Bayesian Ridge Regression Model is now included into our predictiveModeling
#'
#' @export
require(monomvn)
myBridgeModel <- setRefClass(Class = "myBridgeModel",
                             contains="PredictiveModel",
                             fields=c("coefficients","lambda"),
                             methods = list(
                               initialize = function(...){
                                 return(.self)
                               },
                               copy = function() {
                                 
                                 result <- myBridgeModel$new()
                                 result$model <- .self$model
                                 
                                 return(result)
                               },
                               
                               customTrain = function(featureData, responseData, filterData = TRUE, mcmcNum = 1000, burnIn = 200,...){
                                 if(filterData == TRUE){
                                   message("filtering data...")
                                   processedData <- filterPredictiveModelData(featureData, responseData)
                                   featureData <- processedData$featureData
                                   responseData <- processedData$responseData
                                 }
                                 
                                 bridgeModel <- bridge(featureData,responseData,T=mcmcNum,verb=0)
                                 
                                 BURNIN <- burnIn
                                 # discard first 10 ~ 20% of MCMC sampling as a burning-in process
                                 
                                 betaCoefs <- apply(bridgeModel$beta[(BURNIN +1):nrow(bridgeModel$beta),],2,mean)
                                 
                                 # names(betaCoefs) <- colnames(featureData)
                                 
                                 intercept <- mean(bridgeModel$mu[(BURNIN+1):length(bridgeModel$mu)])
                                 betaCoefs <- c(intercept,betaCoefs)
                                 names(betaCoefs)[2:length(betaCoefs)] <- colnames(featureData)
                                 
                                 .self$coefficients <- betaCoefs
                                 .self$lambda <- bridgeModel$lambda2[(BURNIN +1):length(bridgeModel$lambda2)]
                                 
                               },
                               
                               customPredict = function(featureData){
                                 # Be careful when p>>n : This is singular case. So we have to consider "NA" for predicting response
                                 featureNames <- names(.self$coefficients[2:length(.self$coefficients)])
                                 coefNames <-names(.self$coefficients[which(!is.na(.self$coefficients))])
                                 interNames <- intersect(featureNames,coefNames)
                                 predictedResponse <- .self$coefficients[1] + featureData[,interNames] %*% .self$coefficients[interNames]
                                 return(predictedResponse)
                               }
                               )
                             )
