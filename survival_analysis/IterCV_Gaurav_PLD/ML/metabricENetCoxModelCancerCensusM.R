metabricENetCoxModelCancerCensusM <- setRefClass(Class = "metabricENetCoxModelCancerCensusM",                              
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
                                                     
                                                     load("~/COMPBIO/trunk/users/pandey/controlledExptGeneLists.Rdata")
                                                     
                                                     featureData <-createAggregateFeatureDataSet(list(expr=exprData,copy = copyData))
                                                     featureData<-featureData[cancer_census_names,]
                                                     featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                                     featureData <- unique(featureData_filtered)
                                                     
                                                     FEA <-t(featureData)                                                
                                                     
                                                     # Model training
                                                     .self$childclass <- myCoxModel$new()
                                                     .self$model <- .self$childclass$customTrain(FEA,
                                                                                                 clinicalSurvData,
                                                                                                 ...)                                              
                                                     
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
