# alpha = n x 1 vector, lambda = m x 1 vector
# (alpha,lambda) = n x m matrix tune grid
# each alpha elements in alpha vector have cv errors corresponding to lambda vector
require(survival)
crossValidationParameterOptimizer_penalty <- setRefClass(Class = "crossValidationParameterOptimizer_penalty",
                                                         fields=c("cv.all","opt.metricGrid","opt.param","opt.coeffs","model"),
                                                         methods = list(
                                                           
                                                           findOpt = function(model, featureData, responseData, numFolds, alpha, lambda, ...){
                                                             numAlpha = length(alpha)                                                             
                                                             numLambda = length(lambda)                                                             
                                                             
                                                             folds <- createFolds(featureData[,1], k=numFolds)
                             
                                                             alphaModel <- foreach(j = 1:length(alpha)) %dopar% {  
                                                               foldModel <- foreach (k = 1:numFolds) %dopar% {
                                                                 curModel <- model$new()
                                                                 curModel$customTrain(featureData[-folds[[k]],], responseData[-folds[[k]]], alpha = alpha[j], lambda = lambda, ...)
                                                                 return(curModel)
                                                               }
                                                               
                                                               foldTestPredictions <- foreach(i=1:numFolds) %do% foldModel[[i]]$customPredict(featureData[folds[[i]],])
                                                               foldTrainPredictions <- foreach(i=1:numFolds) %do% foldModel[[i]]$customPredict(featureData[-folds[[i]],])                                                           
                                                               foldTestObservations <- foreach(i=1:numFolds) %do% responseData[folds[[i]]]                                                               
                                                               foldTrainObservations <- foreach(i=1:numFolds) %do% responseData[-folds[[i]]]
                                                                                                                            
                                                               return(list(foldTestPredictions = foldTestPredictions, foldTrainPredictions = foldTrainPredictions,
                                                                           foldTestObservations = foldTestObservations, foldTrainObservations = foldTrainObservations))
                                                             }                                                             
                                                             .self$cv.all <-alphaModel
                                                             
                                                             
                                                             # by default, MSE is implemented
                                                             trainErrorFunction <- function(obsr,pred){
                                                               return(mean((pred - obsr)^2))
                                                             }
                                                             # if responseData is survival object, then concordance index is used
                                                             # however, because the more the better, I set 1-cindex to find optimal parameter in order to use min function
                                                             if(is.Surv(responseData)){   
                                                               trainErrorFunction <- function(obsr,pred){
                                                                 Index <- concordance.index(pred,obsr[,1],obsr[,2])                                                                 
                                                                 return(1-Index$c.index)
                                                                 #Index <- D.index(pred,obsr[,1],obsr[,2])                                                                 
                                                                 #return(10-Index$d.index)
                                                               }
                                                             }
                                                             
                                                             Errors<-foreach(j = 1:numAlpha) %dopar% {
                                                               testPredictions <- foreach(k = 1:numFolds) %dopar% {
                                                                 metric<-c()
                                                                 for(i in 1:ncol(alphaModel[[j]]$foldTestPredictions[[k]])){
                                                                   metric<-c(metric,trainErrorFunction(alphaModel[[j]]$foldTestObservations[[k]],alphaModel[[j]]$foldTestPredictions[[k]][,i]))
                                                                 }
                                                                 return(metric)
                                                               }
                                                               metricTestPredictions <- do.call("rbind", testPredictions)
                                                               
                                                               trainPredictions <- foreach(k = 1:numFolds) %dopar% {
                                                                 metric<-c()
                                                                 for(i in 1:ncol(alphaModel[[j]]$foldTrainPredictions[[k]])){
                                                                   metric<-c(metric,trainErrorFunction(alphaModel[[j]]$foldTrainObservations[[k]],alphaModel[[j]]$foldTrainPredictions[[k]][,i]))
                                                                 }
                                                                 return(metric)
                                                               }
                                                               metricTrainPredictions <- do.call("rbind", trainPredictions)
                                                               
                                                               return(list(metricTestPredictions=metricTestPredictions,metricTrainPredictions=metricTrainPredictions))                                                               
                                                             }
                                                             
                                                             # Rather than mean, here I use Quartile mean from 1st Quartile and 3rd Quartile
                                                             IQM <- function(vec,na.rm = TRUE){
                                                               a<-quantile(vec,na.rm = TRUE)
                                                               b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
                                                               return(mean(vec[b],na.rm = TRUE))
                                                             }
                                                             
                                                             meanError<-c()                                                             
                                                             for(j in 1:numAlpha){
                                                               meanError<-rbind(meanError,apply(Errors[[j]]$metricTestPredictions,2,IQM))
                                                               #meanError<-rbind(meanError,apply(Errors[[j]]$metricTestPredictions,2,mean))
                                                             }
                                                             .self$opt.metricGrid<-meanError
                                                             
                                                             if(numAlpha ==1)
                                                               alpha.opt = alpha
                                                             else
                                                               alpha.opt = alpha[which.min(apply(meanError,1,min))]
                                                             
                                                             if(numLambda == 1)
                                                               lambda.opt = lambda
                                                             else
                                                               lambda.opt = lambda[which.min(apply(meanError,2,min))]
                                                             
                                                             .self$opt.param = list(alpha.opt = alpha.opt,lambda.opt = lambda.opt)
                                                             
                                                             # process to find coefficients
                                                             bestModel <- model$new()
                                                             bestModel$customTrain(featureData, responseData, alpha = alpha.opt, lambda = lambda.opt, ...)
                                                             .self$opt.coeffs <- bestModel$getCoefficients()
                                                             .self$model <- bestModel
                                                           }
                                                           
                                                           )
                                                         )