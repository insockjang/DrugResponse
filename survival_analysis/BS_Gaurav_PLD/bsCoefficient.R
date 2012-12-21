bsCoefficient <- setRefClass(Class = "bsCoefficient",                              
#                                           contains="PredictiveModel",
                                          fields=c("model","coefBS"),
                                          methods = list(
                                            initialize = function(...){
                                              return(.self)
                                            },
                                            
                                            rawModel = function(){
                                              return(.self$model)
                                            },
                                            customCoef = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, priorName, ...){
                                              load("~/COMPBIO/trunk/users/pandey/controlledExptGeneLists.Rdata")
                                              
                                              featureData <-createAggregateFeatureDataSet(list(expr=exprData,copy = copyData))
                                              if(is.element("cancer_census",priorName)){
                                                featureData<-featureData[cancer_census_names,]                                                       
                                              }
                                              if(is.element("marginal",priorName)){
                                                featureData<-featureData[marginal_association_names,]                                                       
                                              }
                                              if(is.element("metabric",priorName)){
                                                featureData<-featureData[metabric_clustering_names,]                                                       
                                              }
                                              if(is.element("topvaringhiggins",priorName)){
                                                featureData<-featureData[topvarying_higgins_names,]                                                       
                                              }
                                              if(is.element("topvarings",priorName)){
                                                featureData<-featureData[topvarying_names,]                                                       
                                              }
                                              featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                              featureData <- unique(featureData_filtered)
                                              
                                              FEA <-t(featureData)                                                
                                              
                                              # Model training
                                              .self$model <- myCoxModelBS$new()
                                              .self$model$customTrain(FEA,clinicalSurvData,...)
                                              .self$coefBS<-.self$model$customPredict()
                                              
                                            }
                                            )
                                          )
