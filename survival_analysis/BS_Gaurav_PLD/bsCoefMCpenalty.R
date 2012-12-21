bsCoefMCpenalty <- setRefClass(Class = "bsCoefMCpenalty",                              
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
                           
                           FEA1 <-as.data.frame(t(featureData))
                           FEA <- as.matrix(cbind(FEA1,RES1))                                            
                           
                           # Model training
                           .self$model <- myCoxModelBS$new()
                           .self$model$customTrain(FEA,clinicalSurvData,...)
                           .self$coefBS<-.self$model$customPredict()
                           
                         }
                         )
                       )
