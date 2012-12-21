require(randomForest)
require(modeest)

myRandomForestModel <- setRefClass(Class = "myRandomForestModel",
                        contains="PredictiveModel",
                        fields="model",
                        methods = list(
                          initialize = function(...){
                            return(.self)
                          },
                          
                          rawModel = function(){
                            return(.self$model)
                          },
                          
                          customTrain = function(featureData, responseData,ntree = 50,...){
                            if(!is.factor(responseData)){                              
                              responseData<-factor(as.numeric(responseData<=mlv(responseData)$M))                              
                            }                            
                            .self$model <- randomForest(featureData, responseData, ntree = ntree, do.trace = 5,...)
                          },
                          
                          customPredict = function(featureData){
                            predictedResponse <- predict(.self$model, featureData, type = "prob")
                            return(predictedResponse[,2])
                          }
                          )
                        )
