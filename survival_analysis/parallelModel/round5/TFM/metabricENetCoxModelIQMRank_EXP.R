source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")

metabricENetCoxModelIQMRank_EXP <- setRefClass(Class = "metabricENetCoxModelIQMRank_EXP",                              
                                           fields=c("model","childclass"),
                                           methods = list(
                                             initialize = function(...){
                                               return(.self)
                                             },
                                             
                                             rawModel = function(){
                                               return(.self$model)
                                             },
                                             
                                             customTrain = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, ...){
                                               
                                               load("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TFM_info/TFM.Rdata")
                                               TFM1 <- exprTFM$TF|exprTFM$M
                                               TFM2 <- copyTFM$TF|copyTFM$M
                                               
                                               featureData <-createAggregateFeatureDataSet(list(expr=exprData[TFM1,]))
                                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                               featureData <- unique(featureData_filtered)
                                               
                                               FEA <-t(featureData)                                                
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               # Model training
                                               .self$childclass <- myEnetCoxModelIQM$new()
                                               .self$model <- .self$childclass$customTrain(FEA,
                                                                                           clinicalSurvData,
                                                                                           alpha = alphas, 
                                                                                           lambda = lambdas,
                                                                                           nfolds =10)
                                             },
                                             customPredict = function(exprData, copyData, clinicalFeaturesData, ...){
                                               beta <- rownames(.self$childclass$getCoefficients())
                                               
                                               featureData <-createAggregateFeatureDataSet(list(expr=exprData))
                                               featureData <- t(featureData)
                                               FEA <- as.data.frame(featureData)
                                               FEA <-as.matrix(FEA[,beta])
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               
                                               predictedResponse <- predict(.self$childclass$model,FEA)
                                               return(predictedResponse)
                                             }
                                             )
                                           )

