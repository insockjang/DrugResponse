source("~/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")

M_ExpCNV <- setRefClass(Class = "M_ExpCNV",                              
                        contains="PredictiveModel",
                        fields=c("model","childclass"),
                        methods = list(
                          initialize = function(...){
                            return(.self)
                          },
                          
                          rawModel = function(){
                            return(.self$model)
                          },
                          
                          customTrain = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, ...){
                            
                            featureData <-createAggregateFeatureDataSet(list(expr=exprData,copy = copyData))
                            
                            featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                            featureData <- unique(featureData_filtered)
                            
                            FEA <-t(featureData)                                                
                            
                            # Model training
                            .self$childclass <- myEnetCoxModel$new()
                            .self$model <- .self$childclass$customTrain(FEA,
                                                                        clinicalSurvData,
                                                                        alpha = alphas, 
                                                                        lambda = lambdas,
                                                                        nfolds =10)                                                
                            
                          },
                          customPredict = function(exprData, copyData, clinicalFeaturesData, ...){
                            beta <- rownames(.self$childclass$getCoefficients())
                            featureData <-createAggregateFeatureDataSet(list(expr=exprData,copy = copyData))
                            
                            
                            FEA<-featureData[beta,]
                            FEA <- t(FEA)                                                                                    
                            
                            predictedResponse <- predict(.self$childclass$model,FEA)
                            return(predictedResponse)
                          }
                          )
                        )
