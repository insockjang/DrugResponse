source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")

metabricENetCoxModelMeanRank <- setRefClass(Class = "metabricENetCoxModelMeanRank",                              
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
                                               
                                               featureData <-createAggregateFeatureDataSet(list(expr=exprData[exprTFM$TF,],copy = copyData[copyTFM$TF,]))
                                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                               featureData <- unique(featureData_filtered)
                                               
                                               FEA <-t(featureData)
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               
                                               #bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                               #alphas = seq(0.1,1,len = 10)
                                               #lambdas  = seq(0.001,0.2,len =100)
                                               
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
                                               featureData <- t(featureData)
                                               FEA <- as.data.frame(featureData)
                                               FEA <-as.matrix(FEA[,beta])
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               
                                               predictedResponse <- predict(.self$childclass$model,FEA)
                                               return(predictedResponse)
                                             }
                                             )
                                           )

