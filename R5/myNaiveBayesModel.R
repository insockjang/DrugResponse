require(e1071)
myNaiveBayesModel <- setRefClass(Class = "myNaiveBayesModel",
                                   contains="PredictiveModel",
                                   fields="model",
                                   methods = list(
                                     initialize = function(...){
                                       return(.self)
                                     },
                                     
                                     rawModel = function(){
                                       return(.self$model)
                                     },
                                     
                                     customTrain = function(featureData, responseData,nfolds = NULL){
                                       if(!is.factor(responseData)){
                                         error("response must be binary factor")
                                         break
                                       }                            
                                       .self$model <- naiveBayes(featureData,responseData)
                                     },
                                     
                                     customPredict = function(featureData){
                                       predictedResponse <- predict(.self$model, featureData)
                                       return(predictedResponse)
                                     }
                                     )
                                   )
