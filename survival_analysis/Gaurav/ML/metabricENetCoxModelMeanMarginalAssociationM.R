source("~/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")

metabricENetCoxModelMeanMarginalAssociationM <- setRefClass(Class = "metabricENetCoxModelMeanMarginalAssociationM",                              
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
                                                     featureData<-featureData[marginal_association_names,]
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
