require(survival)
bestCCModel <- setRefClass(Class    = "bestCCModel",
                                 contains = "PredictiveModel",
                                 fields   = c("model"),
                                 methods  = list(
                                   
                                   initialize = function(){                                                                   
                                     return(.self)
                                   },
                                   
                                   copy = function() {
                                     result <- bestCCModel$new()					
                                     result$model  <- .self$model
                                     return(result)
                                   },
                                   
                                   customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...)
                                   {                                     
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
                                     
                                     bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                     .self$model <- coxph(clinicalSurvData ~ .,RES1[,which(bestCC ==1)])
                                     
                                     
                                   },
                                   
                                   customPredict = function(exprData, copyData, clinicalFeaturesData){
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
                                     
                                     bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                     
                                     predict(.self$model, RES1[,which(bestCC ==1)])
                                   }
                                   )
                                 )

