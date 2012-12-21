#' Perform cross validation on a preditive model, this model support Cox model and predictiveModeling
#'
#' @param featureData
#' @param responseData
#' @param model either an instance of PredictiveModel or a string holding the name of one of the machine learning methods that caret supports
#' @param numFolds defaults to 5
#' @param trControl defaults to defaultTrainControl
#' @return a list of PredictiveModelPerformance one per fold
#' @seealso defaultTrainControl
#' @export
findOptimalCoxParam_DIndex <- function(featureData, responseData, model, Grid = createENetTuneGrid(alphas=1),numFolds = 5,...){
  
  #-----------------------------------------------------------------------
  # Make Grid to find optimal alpha and lambda with concordance index as its metric
  # -----------------------------------------------------------------------
  #foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
  
  Coefs <- foreach(grids = 1:nrow(Grid)) %dopar% {
    gridModel <- model$copy()      
    cvResults_cox <- crossValidateFindOptimalCoxModel(featureData, responseData, model = gridModel, numFolds = numFolds, GRID = Grid[grids,])
    #gridModel$customTrain(featureData, responseData, alpha = Grid[grids,1],lambda = Grid[grids,2])
    dindexTrain<-c()
    dindexTest<-c()
    for (i in 1:length(cvResults_cox[[1]])){
      dTrain<-D.index(x=cvResults_cox$trainPredictions[[i]], surv.time=cvResults_cox$trainObservations[[i]][,1], surv.event=cvResults_cox$trainObservations[[i]][,2], na.rm=TRUE, alpha= .05)
      dindexTrain <- c(dindexTrain,dTrain$d.index)
      dTest<-D.index(x=cvResults_cox$testPredictions[[i]], surv.time=cvResults_cox$testObservations[[i]][,1], surv.event=cvResults_cox$testObservations[[i]][,2], na.rm=TRUE, alpha= .05)
      dindexTest <- c(dindexTest,dTest$d.index)
    }
    return(list(CITrain = dindexTrain, CITest=dindexTest))            
  }
  return(Coefs)
}
