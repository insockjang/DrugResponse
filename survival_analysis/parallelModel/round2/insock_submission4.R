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



###################################################
### Call my predictive model class
###################################################
lassoMFbestCCModel <- setRefClass(Class    = "lassoMFbestCCModel", 
                                    contains = "PredictiveModel",
                           fields   = c("model"),
                           methods  = list(
                             
                             initialize = function(){                                                                   
                               return(.self)
                             },
                                                          
                             customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData, ...)
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
                               
                               RES1<-as.matrix(res)             
                               
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               
                               FEA <- cbind(featureData,RES1)
                                                              
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               
                               lambdas = exp(seq(-5,0,len =10))
                               myTest<-myEnetCoxModel$new()
                               .self$model <- myTest$customTrain(FEA,clinicalSurvData,alpha = 1,lambda=lambdas[-c(1:4)],nfolds =5,penalty.factor = c(rep(1,ncol(featureData)),1-bestCC))

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
                               
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               FEA <-cbind(featureData,RES1)
                               
                               return(predict(.self$model,FEA))
                             }
                             
                             )
                           )

###################################################
### trainModel
###################################################

JISModel<- lassoMFbestCCModel$new()
JISModel$customTrain(exprData, copyData, clinicalFeaturesData,clinicalSurvData)
trainPredictions <- JISModel$customPredict(exprData, copyData, clinicalFeaturesData)

###################################################
### submitModel
###################################################
submittedModelParentId <- "161036"
myModelName <- "In Sock Lasso MF plus best CC model"
submittedModelEntity <- Layer(list(name=myModelName, type="E", parentId=submittedModelParentId))
submittedModelEntity <- addObject(submittedModelEntity, JISModel)
submittedModelEntity <- storeEntity(submittedModelEntity)
