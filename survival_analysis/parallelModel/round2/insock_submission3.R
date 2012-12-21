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
models<-loadEntity(162658)
myModel1<-models$objects$myModel1
BETA <- myModel1$getCoefficients()
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
mylassoMFbestCCModel <- setRefClass(Class    = "mylassoMFbestCCModel", 
                                    contains = "PredictiveModel",
                           fields   = c("model"),
                           methods  = list(
                             
                             initialize = function(){                                                                   
                               return(.self)
                             },
                                                          
                             customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,betanames, ...)
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
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               lambda = exp(-1.6666667)
                               dataSets_MF_clinical <- createFeatureAndResponseDataList(featureData, clinicalSurvData)
                               FEA<-as.matrix(cbind(dataSet_MF_clinical$featureData,RES1))
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               .self$model <- glmnet(FEA,clinicalSurvData,family = "cox",alpha =1,lambda
                                                    as.data.frame(cbind(dataSets_MF_clinical$featureData[,betanames],RES1[,which(bestCC!=0)])))
                               
                             },
                             
                             customPredict = function(exprData, copyData, clinicalFeaturesData, betanames,...){
                               
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
                               
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               
                               dataSets_MF_clinical <- createFeatureAndResponseDataList(featureData, clinicalFeaturesData)
                               predict(.self$model,as.data.frame(cbind(dataSets_MF_clinical$featureData[,betanames],RES1[,which(bestCC !=0)])))
                             }                             
                             )
                           )

###################################################
### trainModel
###################################################

JISModel<- mylassoMFbestCCModel$new()
JISModel$customTrain(exprData, copyData, clinicalFeaturesData,clinicalSurvData,betanames = betanames)
trainPredictions <- JISModel$customPredict(exprData, copyData, clinicalFeaturesData,betanames = betanames)


###################################################
### computeTrainCIndex
###################################################
cIndex_train <- concordance.index(x=trainPredictions, surv.time=clinicalSurvData[,"time"], surv.event=clinicalSurvData[,"status"], na.rm=TRUE, alpha= .05)
cIndex_train$c.index # returns the cindex1
cIndex_train$lower # lower CI bound of cindex1
cIndex_train$upper # upper CI bound of cindex1

###################################################
### submitModel
###################################################
submittedModelParentId <- "162653"
myModelName <- "Lasso penalized MF and unpenalized best CC model"

analysis <- Dataset(list(
  name= "Lasso penalized MF and unpenalized best CC model with 5 fold CV",
  description="metric is concordance index when cv is applied",
  parentId="162653"
  ))
(dataset <- createEntity(analysis))

datasetID<-propertyValue(dataset, "id")



submittedModelEntity <- Layer(list(name="pMFuCC", type="E", parentId=datasetID))
submittedModelEntity <- addObject(submittedModelEntity, myTest)
submittedModelEntity <- storeEntity(submittedModelEntity)
