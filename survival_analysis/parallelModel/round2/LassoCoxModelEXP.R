LassoCoxModelEXP <- setRefClass(Class = "LassoCoxModelEXP",
                                contains="PredictiveModel",
                                fields="model",
                                methods = list(
                                  initialize = function(...){
                                    return(.self)
                                  },
                                  
                                  rawModel = function(){
                                    return(.self$model)
                                  },
                                  
                                  customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...){
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
                                    
                                    featureData <-createAggregateFeatureDataSet(list(expr = exprData))
                                    featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                    featureData <- unique(featureData_filtered)
                                    featureData <- scale(t(featureData))
                                    
                                    FEA1 <-as.data.frame(featureData)
                                    FEA <- as.matrix(cbind(FEA1,RES1))
                                    
                                    foldIndices = createFolds(featureData[,1],k = nfolds)
                                    
                                    results<-c()
                                    for(fold in foldIndices){
                                      fit<-glmnet(as.matrix(FEA)[-fold,],clinicalSurvData[-fold],family = "cox", alpha = alpha, lambda = lambdas,penalty.factor = c(rep(1,ncol(featureData)),1-(bestCC)))
                                      pred<-predict(fit,FEA[fold,])
                                      cIndex<-c()
                                      for(k in 1:ncol(pred)){
                                        cIndex<-c(cIndex,concordance.index(pred[,k],clinicalSurvData[fold,1],clinicalSurvData[fold,2])$c.index)
                                      }
                                      results <- rbind(results, cIndex)
                                    }
                                    
                                    
                                    IQM <- function(vec,na.rm = TRUE){
                                      a<-quantile(vec,na.rm = TRUE)
                                      b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
                                      return(mean(vec[b],na.rm = TRUE))
                                    }
                                    
                                    #metrics <- apply(results,2,mean)
                                    metrics <- apply(results,2,IQM)
                                    optParam <- c(max(metrics), alpha, fit$lambda[which.max(metrics)])
                                    names(optParam) <- c("cIndex","alpha","lambda")
                                    .self$model <- glmnet(FEA,clinicalSurvData,family = "cox", alpha = optParam[2], lambda = optParam[3],penalty.factor = c(rep(1,ncol(featureData)),1-(bestCC)))
                                    .self$model$optParam <- optParam 
                                    .self$model$results <- results
                                  },
                                  
                                  customPredict = function(exprData, copyData, clinicalFeaturesData, ...){
                                    beta<-rownames(coef(.self$model))
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
                                    
                                    featureData <-createAggregateFeatureDataSet(list(expr = exprData))
                                    featureData <- scale(t(featureData))
                                    FEA1 <- as.data.frame(featureData)
                                    FEA <-cbind(FEA1,RES1)
                                    FEA <-as.matrix(FEA[,beta])
                                    
                                    predictedResponse <- predict(.self$model,FEA)
                                    return(predictedResponse)
                                  },
                                  getCoefficients = function(){
                                    return(coef(.self$model))                              
                                  }                                                       
                                  
                                  )
                                )

