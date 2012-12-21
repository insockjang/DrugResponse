metabricENetCoxModelMprior <- setRefClass(Class = "metabricENetCoxModelpriorM",                              
                                                 contains="PredictiveModel",
                                                 fields=c("model"),
                                                 methods = list(
                                                   initialize = function(...){
                                                     return(.self)
                                                   },
                                                   
                                                   rawModel = function(){
                                                     return(.self$model)
                                                   },
                                                   
                                                   customTrain = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, priorName, beta,...){
                                                     
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
                                                                                                                                                               
                                                     .self$model<-coxph(clinicalSurvData ~.,data = as.data.frame(FEA[,beta]))  
                                                     .self$model$beta<-beta
                                                     
                                                     
                                                   },
                                                   customPredict = function(exprData, copyData, clinicalFeaturesData,...){
                                                     beta<-.self$model$beta
                                                     featureData <-createAggregateFeatureDataSet(list(expr=exprData,copy = copyData))
                                                                                                          
                                                     FEA<-featureData[beta,]
                                                     FEA <- t(FEA)                                                                                    
                                                     
                                                     predictedResponse <- predict(.self$model,as.data.frame(FEA))
                                                     return(predictedResponse)
                                                   }
                                                   )
                                                 )
