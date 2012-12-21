#' Constructor for a class in the PredictiveModel class hierarchy that implements a simple function regressing the top correlated features
#'
#' @export
myLsModel <- setRefClass(Class = "myLsModel",
                                      contains="PredictiveModel",
                                      fields=list(coefficients="numeric"),
                                      methods = list(
                                        initialize = function(...){
                                          return(.self)
                                        },
                                        
                                        train = function(featureData, responseData, filterData = TRUE, ...){
                                          if(filterData == TRUE){
                                            message("filtering data...")
                                            processedData <- filterPredictiveModelData(featureData, responseData)
                                            featureData <- processedData$featureData
                                            responseData <- processedData$responseData
                                          }
                                          print(dim(featureData))
                                          print(length(responseData))
                                          
                                          
                                          fitCoefs <- lm(responseData ~ (featureData))$coefficients
                                          names(fitCoefs)[2:length(fitCoefs)] <- colnames(featureData)                                          
                                          .self$coefficients <- fitCoefs
                                        },
                                        
                                        predict = function(featureData){
                                          # Be careful when p>>n : This is singular case. So we have to consider "NA" for predicting response
                                          featureNames <- names(.self$coefficients[2:length(.self$coefficients)])
                                          coefNames <-names(.self$coefficients[which(!is.na(.self$coefficients))])
                                          interNames <- intersect(featureNames,coefNames)
                                          predictedResponse <- .self$coefficients[1] + featureData[,interNames] %*% .self$coefficients[interNames]
                                          return(predictedResponse)
                                        }
                                        )
                                      )
