###################################################
### loadLibraries and predefined classes
###################################################
library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(glmnet)
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelMean.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelIQM.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelMeanRank.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelIQMRank.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelMeanRank_CNV.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelIQMRank_CNV.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelMeanRank_EXP.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round2/metabricENetCoxModelIQMRank_EXP.R")


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

bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
alphas = seq(0.1,1,len = 10)
lambdas  = seq(0.001,0.2,len =100)

myMetabricENetCoxModelMean_CNV <- metabricENetCoxModelMean$new()
myMetabricENetCoxModelMean_CNV$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMean_CNV$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)

myMetabricENetCoxModelIQM_CNV <- metabricENetCoxModelIQM$new()
myMetabricENetCoxModelIQM_CNV$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictionsIQM <- myMetabricENetCoxModelIQM_CNV$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)

myMetabricENetCoxModelMean <- metabricENetCoxModelMean$new()
myMetabricENetCoxModelMean$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMean$customPredict(exprData, copyData, clinicalFeaturesData)
testPredictions <- myMetabricENetCoxModelMean$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)

myMetabricENetCoxModelIQM <- metabricENetCoxModelIQM$new()
myMetabricENetCoxModelIQM$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictionsIQM <- myMetabricENetCoxModelIQM$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)

myMetabricENetCoxModelMeanRank <- metabricENetCoxModelMeanRank$new()
myMetabricENetCoxModelMeanRank$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMeanRank$customPredict(exprData, copyData, clinicalFeaturesData)
testPredictions <- myMetabricENetCoxModelMeanRank$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)

myMetabricENetCoxModelIQMRank <- metabricENetCoxModelIQMRank$new()
myMetabricENetCoxModelIQMRank$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelIQMRank$customPredict(exprData, copyData, clinicalFeaturesData)
testPredictions <- myMetabricENetCoxModelIQMRank$customPredict(exprData_validation, copyData_validation, clinicalFeaturesData_validation)


myMetabricENetCoxModelMeanRank_CNV <- metabricENetCoxModelMeanRank_CNV$new()
myMetabricENetCoxModelMeanRank_CNV$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMeanRank_CNV$customPredict(exprData, copyData, clinicalFeaturesData)

myMetabricENetCoxModelIQMRank_CNV <- metabricENetCoxModelIQMRank_CNV$new()
myMetabricENetCoxModelIQMRank_CNV$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelIQMRank_CNV$customPredict(exprData, copyData, clinicalFeaturesData)

myMetabricENetCoxModelMeanRank_EXP <- metabricENetCoxModelMeanRank_EXP$new()
myMetabricENetCoxModelMeanRank_EXP$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMeanRank_EXP$customPredict(exprData, copyData, clinicalFeaturesData)

myMetabricENetCoxModelIQMRank_EXP <- metabricENetCoxModelIQMRank_EXP$new()
myMetabricENetCoxModelIQMRank_EXP$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelIQMRank_EXP$customPredict(exprData, copyData, clinicalFeaturesData)



cIndex_train <- concordance.index(x=trainPredictions, surv.time=clinicalSurvData[,"time"], surv.event=clinicalSurvData[,"status"], na.rm=TRUE, alpha= .05)
cIndex_train$c.index # returns the cindex1
cIndex_train$lower # lower CI bound of cindex1
cIndex_train$upper # upper CI bound of cindex1



submittedModelParentId <- "161036"
myModelName <- "In Sock RANK Mean ENet penalty EXP and unpenalty best CC"
submittedModelEntity <- Layer(list(name=myModelName, type="E", parentId=submittedModelParentId))
submittedModelEntity <- addObject(submittedModelEntity, myMetabricENetCoxModelMeanRank_EXP)
submittedModelEntity <- storeEntity(submittedModelEntity)

