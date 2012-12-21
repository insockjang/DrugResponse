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
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               
                               dataSets_MF_clinical <- createFeatureAndResponseDataList(featureData, clinicalSurvData)
                               dataSets_CC_clinical <- createFeatureAndResponseDataList(pData(clinicalFeaturesData), clinicalSurvData)
                               
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               .self$model <- coxph(clinicalSurvData ~ .,
                                                    as.data.frame(cbind(dataSets_MF_clinical$featureData[,betanames],dataSets_CC_clinical$featureData[,which(bestCC!=0)])))
                               
                             },
                             
                             customPredict = function(exprData, copyData, clinicalFeaturesData, betanames,...){
                               featureData <-createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))
                               featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "rows")
                               featureData <- unique(featureData_filtered)
                               featureData <- scale(t(featureData))
                               
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               
                               dataSets_MF_clinical <- createFeatureAndResponseDataList(featureData, clinicalFeaturesData)
                               predict(.self$model,as.data.frame(cbind(dataSets_MF_clinical$featureData[,betanames],dataSets_MF_clinical$responseData[,which(bestCC !=0)])))
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
### submitModel
###################################################
submittedModelParentId <- "161036"
myModelName <- "In Sock Lasso MF plus best CC model"
submittedModelEntity <- Layer(list(name=myModelName, type="E", parentId=submittedModelParentId))
submittedModelEntity <- addObject(submittedModelEntity, JISModel)
submittedModelEntity <- storeEntity(submittedModelEntity)
