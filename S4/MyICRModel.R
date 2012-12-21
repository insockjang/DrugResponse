#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
MyICRModel <- setRefClass(Class = "MyICRModel",
                          contains="PredictiveModel",
                          fields=c("model", "modelType"),
                          methods = list(
                            initialize = function(modelType, ...){
                              .self$modelType = modelType
                              
                              return(.self)
                            },
                            
                            rawCaretModel = function(){
                              return(.self$model)
                            },
                            
                            train = function(featureData, responseData, trControl = defaultTrainControl(),
                                             filterData = TRUE, tuneGrid = expand.grid(.n.comp = c(1,2,3,4)){
                              if(filterData == TRUE){
                                processedData <- filterPredictiveModelData(featureData, responseData)
                                featureData <- processedData$featureData
                                responseData <- processedData$responseData
                              }
                              .self$model <- caret::train(featureData, 
                                                          responseData, 
                                                          method = .self$modelType,
                                                          preProcess = NULL, # preProcess = c("center", "scale"),
                                                          tuneLength = 4,
                                                          trControl = trControl,
                                                          tuneGrid = tuneGrid)
                                                          ) 
                            },
                            
                            predict = function(featureData){
                              predictedResponse <- predict.train(.self$model, featureData)
                              return(predictedResponse)
                            }
                            )
                          )
