###################################################
### loadLibraries and predefined classes
###################################################
library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(glmnet)

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

source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelExpCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelExp.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelExpCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelCnv.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelExp.R")

myRankLassoExpCnv <- rankLassoModelExpCnv$new()
myRankLassoExpCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoExpCnv$customPredict(exprData, copyData, clinicalFeaturesData)

myRankLassoExp <- rankLassoModelExp$new()
myRankLassoExp$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoExp$customPredict(exprData, copyData, clinicalFeaturesData)


myRankLassoCnv <- rankLassoModelCnv$new()
myRankLassoCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankLassoCnv$customPredict(exprData, copyData, clinicalFeaturesData)

myRank RidgeExpCnv <- rankRidgeModelExpCnv$new()
myRankRidgeExpCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeExpCnv$customPredict(exprData, copyData, clinicalFeaturesData)


myRankRidgeExp <- rankRidgeModelExp$new()
myRankRidgeExp$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeExp$customPredict(exprData, copyData, clinicalFeaturesData)

myRankRidgeCnv <- rankRidgeModelCnv$new()
myRankRidgeCnv$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData)
trainPredictions <- myRankRidgeCnv$customPredict(exprData, copyData, clinicalFeaturesData)


submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Lasso EXP CNV Model", 
                       trainedModel= myRankLassoExpCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelExpCnv.R", 
                       classFileName="rankLassoModelExpCnv")
submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Lasso EXP Model", 
                       trainedModel= myRankLassoExp,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelExp.R", 
                       classFileName="rankLassoModelExp")
submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Lasso CNV Model", 
                       trainedModel= myRankLassoCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankLassoModelCnv.R", 
                       classFileName="rankLassoModelCnv")
submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Ridge EXP CNV Model", 
                       trainedModel= myRankRidgeExpCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelExpCnv.R", 
                       classFileName="rankRidgeModelExpCnv")
submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Ridge EXP Model", 
                       trainedModel= myRankRidgeExp,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelExp.R", 
                       classFileName="rankRidgeModelExp")
submitCompetitionModel(modelName = "In Sock Serial Rank TF1 Ridge CNV Model", 
                       trainedModel= myRankRidgeCnv,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TF/rankRidgeModelCnv.R", 
                       classFileName="rankRidgeModelCnv")
