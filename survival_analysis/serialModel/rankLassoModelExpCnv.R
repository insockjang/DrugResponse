rankLassoModelExpCnv <- setRefClass(Class = "rankLassoModelExpCnv",
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
                                        
                                        featureData <-createAggregateFeatureDataSet(list(expr = exprData,copy = copyData))
                                        featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                        featureData <- unique(featureData_filtered)
                                        
                                        FEA <-as.data.frame(t(featureData))
                                        FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
                                        
                                        alphas = seq(0.1,1,len =10)
                                        lambdas  = exp(seq(-5,0,len =25))
                                        
                                        results <- list()
                                        allFeature <- c()
                                        for(k in 1:14) {
                                          result<-list()
                                          for(k1 in 1 :length(alphas)) {
                                            
                                            exprCnvFit<-cv.glmnet(FEA,RES1[,k],alpha=alphas[k1],lambda = lambdas,nfolds =10)
                                            cvm <- min(exprCnvFit$cvm)
                                            lambda <- exprCnvFit$lambda[which.min(exprCnvFit$cvm)]
                                            
                                            result[[k1]]<-(list(fit =exprCnvFit,cvm = cvm, lambda = lambda))
                                          }
                                          
                                          cv <- c()
                                          for(k1 in 1:length(alphas)){
                                            cv<-c(cv,result[[k1]]$cvm)
                                          }
                                          optimalNum<-which.min(cv)
                                          results[[k]]<-result[[optimalNum]]
                                          Beta<-coef(result[[optimalNum]]$fit,s="lambda.min")
                                          allFeature<-union(allFeature,which(Beta!=0))  
                                        }
                                        
                                        
                                        FEA1<-as.data.frame(FEA)
                                        
                                        .self$model<-cv.glmnet(as.matrix(FEA1[,(allFeature[-1]-1)]),clinicalSurvData,family = "cox",nfolds =10)
                                      },
                                      customPredict=function(exprData,copyData,clinicalFeaturesData,clinicalSurvData){
                                        
                                        beta <- rownames(coef(.self$model,s="lambda.min"))
                                        
                                        featureData <-createAggregateFeatureDataSet(list(expr = exprData,copy = copyData))
                                        featureData <- t(featureData)
                                        FEA1 <- as.data.frame(featureData)
                                        FEA <-as.matrix(FEA1[,beta])
                                        FEA <- as.matrix((apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5)
                                        
                                        
                                        predictedResponse <- predict(.self$model,FEA,s="lambda.min")
                                        return(predictedResponse)                                        
                                      }
                                      )
                                    )
