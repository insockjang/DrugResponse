# genomic instability index with molecular feature only
metabricENetCoxModelGIIM <- setRefClass(Class = "metabricENetCoxModelGIIM",                              
                                        contains="PredictiveModel",
                                        fields=c("model","childclass"),
                                        methods = list(
                                          initialize = function(...){
                                            return(.self)
                                          },
                                          
                                          rawModel = function(){
                                            return(.self$model)
                                          },
                                          
                                          customTrain = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, GII = 1, ...){
                                            
                                            Ginstab = unlist(apply(exprs(copyData), 1, function(x){return(sum(abs(x)>GII, na.rm = TRUE) ) } ) )
                                            GII_name<-paste(names(Ginstab[which(Ginstab>0)]),"_copy",sep="")
                                            
                                            featureData <-createAggregateFeatureDataSet(list(copy = copyData))
                                            featureData<-featureData[GII_name,]
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
                                            featureData <-createAggregateFeatureDataSet(list(copy = copyData))
                                            
                                            
                                            FEA<-featureData[beta,]
                                            FEA <- t(FEA)                                                                                    
                                            
                                            predictedResponse <- predict(.self$childclass$model,FEA)
                                            return(predictedResponse)
                                          }
                                          )
                                        )
