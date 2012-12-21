#' Perform cross validation on a ENet model in order to find optimal alpha and lambdas in predictiveModel packages
#'
#' @param featureData
#' @param responseData
#' @param model R5 class 
#' @param numFolds defaults to 5
#' @param trControl defaults to defaultTrainControl
#' @return a list of PredictiveModelPerformance one per fold
#' @seealso defaultTrainControl
#' @export

cvENetModel <- function(featureData, responseData, numFolds = numFolds, alpha = alpha, lambda = lambda, ...){
  
    
    #-----------------------------------------------------------------------
    # Split the data into training and test partitions
    # -----------------------------------------------------------------------
    
    #set.seed(2)
    foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
    
    message("Training partitions")
    
    modelResults <- foreach(fold = foldIndices) %dopar% {
      foldModel <- ENetModel$new()
      foldModel$customTrain(featureData[-fold,], responseData[-fold], alpha = alpha, lambda = lambda)
      return(foldModel)
    }
    
    message("finding Optimal parameters")
    
    foldTestPredictions <- foreach(i=1:numFolds) %do% modelResults[[i]]$customPredict(featureData[foldIndices[[i]],])
    
    foldTestObservations <- foreach(i=1:numFolds) %do% responseData[foldIndices[[i]]]
    
    RMSE<-c()
    
    for(i in 1:numFolds){
      rmse <-c()
      for(j in 1:dim(foldTestPredictions[[1]])[2]){
        rmse <- c( rmse,mean((foldTestPredictions[[i]][,j]-foldTestObservations[[i]])^2))
      }
      RMSE<-rbind(RMSE,rmse)
    }
    return(RMSE)
}
