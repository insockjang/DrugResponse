#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @export
CrossValidationParameterOptimizer <- setRefClass(Class = "CrossValidationParameterOptimizer",
                                                 contains="PredictiveModel",
                                                 fields=c("cvErrors","optParam"),
                                                 methods = list(
                                                   
                                                   findOptimalParameters = function(model, featureData, responseData, numFolds, tuneGrid,...){
                                                     # if(is(responseData,Surv")
                                                     trainErrorFunction <- function(pred,obsr){
                                                       return(sum((pred - obsr)^2))
                                                     }
                                                     
                                                     folds <- createFolds(featureData[,1], k=numFolds)
                                                     
                                                     cvModel <-foreach(tuneRow = 1:nrow(tuneGrid)) %dopar% {
                                                       
                                                       foldModel <- foreach (fold = folds) %dopar% {
                                                         curModel <- model$copy()
                                                         curModel$customTrain(featureData[-fold,], responseData[-fold], alpha = tuneGrid[tuneRow,1], lambda= tuneGrid[tuneRow,2],...)
                                                         return(curModel)
                                                       }
                                                       foldPredictions<-foreach(i = 1:length(foldModel)) %do% {trainErrorFunction(responseData[folds[[i]]],foldModel[[i]]$customPredict(featureData[folds[[i]],]))}
                                                       foldError<-do.call("c",foldPredictions)
                                                       return(list(cvErrors = foldError,cvfoldModel=foldModel))
                                                     }
                                                     cvPredictions<-foreach(i = 1:length(cvModel)) %do% {mean(cvModel[[i]]$cvErrors)}
                                                     cvError<-do.call("c",cvPredictions)
                                                     
                                                     plot(cvError)
                                                     ##### tunable way to select minimum. Default is minimum. Could also be minimum within standard error.
                                                     .self$optParam<-tuneGrid[which.min(cvError),]                                                     
                                                     .self$cvErrors <- cvError
                                                   }
                                                   )
                                                 )
CrossValidationParameterOptimizer_penalty <- setRefClass(Class = "CrossValidationParameterOptimizer_penalty",
                                                         contains="PredictiveModel",
                                                         fields=c("cvErrors","alpha.opt","lambda.opt"),
                                                         methods = list(
                                                           
                                                           findOptimalParameters = function(model, featureData, responseData, numFolds, alpha, lambda,...){
                                                             # if(family == "cox")
                                                             trainErrorFunction <- function(pred,obsr){
                                                               return(mean((pred - obsr)^2))
                                                             }
                                                             
                                                             folds <- createFolds(featureData[,1], k=numFolds)
                                                             
                                                             # alpha = n x 1 vector, lambda = m x 1 vector
                                                             # (alpha,lambda) = n x m matrix tune grid
                                                             # each alpha elements in alpha vector have cv errors corresponding to lambda vector
                                                             
                                                             alpha = unique(alpha)
                                                             lambda = unique(lambda)
                                                             
                                                             alphaModel <- foreach(j = 1:length(alpha)) %dopar%{
                                                               
                                                               foldModel <- foreach (fold = folds) %dopar% {
                                                                 curModel <- model$new()
                                                                 curModel$customTrain(featureData[-fold,], responseData[-fold], alpha = alpha[j], lambda= lambda,...)
                                                                 return(curModel)
                                                               }
                                                               foldPredictions<-foreach(i = 1:length(foldModel)) %do% {
                                                                 errors<-c()
                                                                 for(k in 1:length(lambda)){
                                                                   errors<-c(errors,trainErrorFunction(responseData[folds[[i]]],foldModel[[i]]$customPredict(featureData[folds[[i]],])[,k]))
                                                                 }
                                                                 return(errors)
                                                               }
                                                               foldError<-do.call("rbind",foldPredictions)
                                                               FoldError<-colMeans(foldError)
                                                               return(FoldError)                                                             
                                                             }
                                                             alphaError <-do.call("rbind",alphaModel)
                                                             rownames(alphaError)<-alpha
                                                             colnames(alphaError)<-lambda
                                                             
                                                             #plot(apply(foldError,2,median))                                                      
                                                             #points(FoldError,col = "red")
                                                             ##### tunable way to select minimum. Default is minimum. Could also be minimum within standard error.
                                                             .self$alpha.opt<-alpha[which.min(apply(alphaError,1,min))]                                                     
                                                             .self$lambda.opt<-lambda[which.min(apply(alphaError,2,min))]                                                                                                                  
                                                             .self$cvErrors <- alphaError
                                                           }
                                                           )
                                                         )