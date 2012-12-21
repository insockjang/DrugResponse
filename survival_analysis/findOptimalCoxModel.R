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
findOptimalCoxModel <- function(featureData, responseData, model, Grid = createENetTuneGrid(alphas=1),numFolds = 5,...){
    
    #-----------------------------------------------------------------------
    # Make Grid to find optimal alpha and lambda with concordance index as its metric
    # -----------------------------------------------------------------------
    #foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
  
    Coefs <- foreach(grids = 1:dim(Grid)[1]) %dopar% {
      gridModel <- model$copy()      
      cvResults_cox <- crossValidateFindOptimalCoxModel(featureData, responseData, model = gridModel, numFolds = numFolds, GRID = Grid[grids,])
      #gridModel$customTrain(featureData, responseData, alpha = Grid[grids,1],lambda = Grid[grids,2])
      cindexTrain<-c()
      cindexTest<-c()
      for (i in 1:length(cvResults_cox[[1]])){
        cTrain<-concordance.index(x=cvResults_cox$trainPredictions[[i]], surv.time=cvResults_cox$trainObservations[[i]][,1], surv.event=cvResults_cox$trainObservations[[i]][,2], na.rm=TRUE, alpha= .05,method="noether")
        cindexTrain <- c(cindexTrain,cTrain$c.index)
        cTest<-concordance.index(x=cvResults_cox$testPredictions[[i]], surv.time=cvResults_cox$testObservations[[i]][,1], surv.event=cvResults_cox$testObservations[[i]][,2], na.rm=TRUE, alpha= .05,method="noether")
        cindexTest <- c(cindexTest,cTest$c.index)
      }
      return(list(CITrain = cindexTrain, CITest=cindexTest))            
    }
    return(Coefs)
}
