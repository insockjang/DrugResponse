###################################################
### loadLibraries and predefined classes
###################################################
library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(glmnet)
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelExpCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelExp.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelExpCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelExp.R")

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
#preprocessed results from lasso in order to select Molecular features

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


myRankLassoExpCnv <- rankLassoModelExpCnv$new()
myRankLassoExpCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoExpCnv$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Lasso EXP CNV Model", 
                       trainedModel= myRankLassoExpCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelExpCnv.R", 
                       classFileName="rankLassoModelExpCnv")

myRankLassoExp <- rankLassoModelExp$new()
myRankLassoExp$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoExp$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Lasso EXP Model", 
                       trainedModel= myRankLassoExp,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelExp.R", 
                       classFileName="rankLassoModelExp")

                       
myRankLassoCnv <- rankLassoModelCnv$new()
myRankLassoCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoCnv$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Lasso CNV Model", 
                       trainedModel= myRankLassoCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankLassoModelCnv.R", 
                       classFileName="rankLassoModelCnv")


myRank RidgeExpCnv <- rankRidgModeleExpCnv$new()
myRankRidgeExpCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeExpCnv$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Ridge EXP CNV Model", 
                       trainedModel= myRankRidgeExpCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelExpCnv.R", 
                       classFileName="rankRidgeModelExpCnv")


myRankRidgeExp <- rankRidgeModelExp$new()
myRankRidgeExp$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeExp$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Ridge EXP Model", 
                       trainedModel= myRankRidgeExp,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelExp.R", 
                       classFileName="rankRidgeModelExp")
                       
myRankRidgeCnv <- rankRidgeModelCnv$new()
myRankRidgeCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeCnv$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)
submitCompetitionModel(modelName = "In Sock Serial Rank Ridge CNV Model", 
                       trainedModel= myRankRidgeCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/rankRidgeModelCnv.R", 
                       classFileName="rankRidgeModelCnv")

