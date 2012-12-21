###################################################
### loadLibraries and predefined classes
###################################################
library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(glmnet)
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")

#preprocessed results from lasso in order to select Molecular features
models<-loadEntity(162982)
myTest<-models$objects$myTest
BETA <- myTest$getCoefficients()
betanames<-rownames(BETA)[which(BETA !=0)]

###################################################
### loadData
###################################################
# synapseLogin() ### not required if configured for automatic login

idExpressionLayer <- "160776" ##"160644" 
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idCopyLayer <- "160778" ##"160646"
copyLayer <- loadEntity(idCopyLayer)
copyData <- copyLayer$objects[[1]]

idClinicalFeaturesLayer <- "160780" ##"139171"
clinicalFeaturesLayer <- loadEntity(idClinicalFeaturesLayer)
clinicalFeaturesData <- clinicalFeaturesLayer$objects[[1]]

idClinicalSurvLayer <- "160782"
clinicalSurvLayer <- loadEntity(idClinicalSurvLayer)
clinicalSurvData <- clinicalSurvLayer$objects[[1]]


myLassoCoxModel <- setRefClass(Class = "myLassoCoxModel",
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
                                  
                                  RES1<-as.matrix(res)             
                                  
                                  featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                                  featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                  featureData <- unique(featureData_filtered)
                                  featureData <- scale(t(featureData))
                                  
                                  FEA <- cbind(featureData,RES1)
                                  
                                  bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                  alpha =1
                                  #lambdas = exp(seq(-2,-1,len = 10))
                                  lambdas  = seq(0.001,0.2,len =1000)
                                  nfolds =20
                                  foldIndices = createFolds(featureData[,1],k = nfolds)
                                  
                                  results<-c()
                                  for(fold in foldIndices){
                                    fit<-glmnet(FEA[-fold,],clinicalSurvData[-fold],family = "cox", alpha = alpha, lambda = lambdas,penalty.factor = c(rep(1,ncol(featureData)),1-(bestCC)),maxit=10000000)
                                    pred<-predict(fit,FEA[fold,])
                                    cIndex<-c()
                                    for(k in 1:ncol(pred)){
                                      cIndex<-c(cIndex,survConcordance(clinicalSurvData[fold]~pred[,k])$concordance)
                                    }
                                    results <- rbind(results, cIndex)
                                  }
                                  metrics <- apply(results,2,mean)
                                  optParam <- c(max(metrics), alpha, lambdas[which.max(metrics)])
                                  names(optParam) <- c("cIndex","alpha","lambda")
                                  .self$model <- glmnet(FEA,clinicalSurvData,family = "cox", alpha = optParam[2], lambda = optParam[3],penalty.factor = c(rep(1,ncol(featureData)),1-(bestCC)),maxit=10000000)
                                  .self$model$optParam <- optParam                                                                            
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
                                  
                                  RES1<-as.matrix(res) 
                                  bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                                  RES1<-RES1[,which(bestCC !=0)]
                                 
                                  featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                                  featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                                  featureData <- unique(featureData_filtered)
                                  featureData <- scale(t(featureData))
                                  
                                  FEA <-cbind(featureData,RES1)
                                  
                                  predictedResponse <- predict(.self$model,FEA)
                                  return(predictedResponse)
                                },
                                getCoefficients = function(){
                                  return(coef(.self$model))                              
                                }                                                       
                                
                                )
                              )

###################################################
### trainModel
###################################################

JISModel1<- myLassoCoxModel$new()
JISModel1$customTrain(exprData, copyData, clinicalFeaturesData,clinicalSurvData)
trainPredictions <- JISModel1$customPredict(exprData, copyData, clinicalFeaturesData)


cIndex_train <- concordance.index(x=trainPredictions, surv.time=clinicalSurvData[,"time"], surv.event=clinicalSurvData[,"status"], na.rm=TRUE, alpha= .05)
cIndex_train$c.index # returns the cindex1
cIndex_train$lower # lower CI bound of cindex1
cIndex_train$upper # upper CI bound of cindex1




###################################################
### submitModel
###################################################
submittedModelParentId <- "161036"
myModelName <- "In Sock LASSO penalized MF plus unpenalized best CC model version 06"
submittedModelEntity <- Layer(list(name=myModelName, type="E", parentId=submittedModelParentId))
submittedModelEntity <- addObject(submittedModelEntity, JISModel1)
submittedModelEntity <- storeEntity(submittedModelEntity)