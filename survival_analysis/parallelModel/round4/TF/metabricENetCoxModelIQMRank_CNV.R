
metabricENetCoxModelIQMRank_CNV <- setRefClass(Class = "metabricENetCoxModelIQMRank_CNV",                              
                                           fields=c("model","childclass"),
                                           methods = list(
                                             initialize = function(...){
                                               return(.self)
                                             },
                                             
                                             rawModel = function(){
                                               return(.self$model)
                                             },
                                             
                                             customTrain = function(exprData,copyData,clinicalFeaturesData,clinicalSurvData, ...){
                                               
                                               RES<-pData(clinicalFeaturesData)
                                               # clinical covariate numeric mapping
                                               res<-matrix(0,ncol = ncol(RES),nrow = nrow(RES))
                                               rownames(res)<-rownames(RES)
                                               colnames(res)<-colnames(RES)[1:14]
                                               
                                               res[,"Site"] <- as.numeric(RES[,"Site"])
                                               res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
                                               res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
                                               res[,"grade"] <- as.numeric(RES[,"grade"])
                                               res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                                                                           ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                                                                  ifelse(RES[,"histology"]== "Medullary",3,
                                                                                         ifelse(RES[,"histology"]== "MixedHistology",4,5))))
                                               
                                               res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
                                               res[,"chemo"] <- as.numeric(RES[,"chemo"])
                                               res[,"hormone"] <- as.numeric(RES[,"hormone"])
                                               res[,"radiation"] <- as.numeric(RES[,"radiation"])
                                               res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
                                               res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
                                               res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
                                               res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
                                               res[,"tripleNegative"] <- as.numeric(RES[,"tripleNegative"])
                                               
                                               RES1<-as.data.frame(res)             
                                               load("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TFM_info/TFM.Rdata")
                                               
                                               featureData <-createAggregateFeatureDataSet(list(copy = copyData[copyTFM$TF,]))
                                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                               featureData <- unique(featureData_filtered)
                                               #featureData <- scale(t(featureData))
                                               
                                               #featureData <- featureData[,1:100]
                                               
                                               FEA1 <-as.data.frame(t(featureData))
                                               FEA <- as.matrix(cbind(FEA1,RES1))
                                               
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               
                                               #bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                               #alphas = seq(0.1,1,len = 10)
                                               #lambdas  = seq(0.001,0.2,len =100)
                                               
                                               # Model training
                                               .self$childclass <- myEnetCoxModelIQM$new()
                                               .self$model <- .self$childclass$customTrain(FEA,
                                                                                           clinicalSurvData,
                                                                                           alpha = alphas, 
                                                                                           lambda = lambdas,
                                                                                           nfolds =10)
                                             },
                                             customPredict = function(exprData, copyData, clinicalFeaturesData, ...){
                                               RES<-pData(clinicalFeaturesData)
                                               # clinical covariate numeric mapping
                                               res<-matrix(0,ncol = ncol(RES),nrow = nrow(RES))
                                               rownames(res)<-rownames(RES)
                                               colnames(res)<-colnames(RES)[1:14]
                                               
                                               res[,"Site"] <- as.numeric(RES[,"Site"])
                                               res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
                                               res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
                                               res[,"grade"] <- as.numeric(RES[,"grade"])
                                               res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                                                                           ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                                                                  ifelse(RES[,"histology"]== "Medullary",3,
                                                                                         ifelse(RES[,"histology"]== "MixedHistology",4,5))))
                                               
                                               res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
                                               res[,"chemo"] <- as.numeric(RES[,"chemo"])
                                               res[,"hormone"] <- as.numeric(RES[,"hormone"])
                                               res[,"radiation"] <- as.numeric(RES[,"radiation"])
                                               res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
                                               res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
                                               res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
                                               res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
                                               res[,"tripleNegative"] <- as.numeric(RES[,"tripleNegative"])
                                               
                                               RES1<-as.data.frame(res) 
                                               
                                               beta <- rownames(.self$childclass$getCoefficients())
                                               
                                               featureData <-createAggregateFeatureDataSet(list(copy = copyData))
                                               featureData <- t(featureData)
                                               FEA1 <- as.data.frame(featureData)
                                               FEA <-cbind(FEA1,RES1)
                                               FEA <-as.matrix(FEA[,beta])
                                               FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                               
                                               
                                               predictedResponse <- predict(.self$childclass$model,FEA)
                                               return(predictedResponse)
                                             }
                                             )
                                           )

