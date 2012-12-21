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

source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModelIQM.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank_CNV.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank_CNV.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank_EXP.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank_EXP.R")

bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
alphas = seq(0.1,1,len = 10)
lambdas  = seq(0.001,0.2,len =100)

myMetabricENetCoxModelMeanRank <- metabricENetCoxModelMeanRank$new()
myMetabricENetCoxModelMeanRank$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelMeanRank$customPredict(exprData, copyData, clinicalFeaturesData)

myMetabricENetCoxModelIQMRank <- metabricENetCoxModelIQMRank$new()
myMetabricENetCoxModelIQMRank$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas,bestCC = bestCC)
trainPredictions <- myMetabricENetCoxModelIQMRank$customPredict(exprData, copyData, clinicalFeaturesData)

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



submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet Mean EXP CNV Model", 
                       trainedModel= myMetabricENetCoxModelMeanRank,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank.R", 
                       classFileName="metabricENetCoxModelMeanRank")
submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet Mean CNV Model", 
                       trainedModel= myMetabricENetCoxModelMeanRank_CNV,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank_CNV.R", 
                       classFileName="metabricENetCoxModelMeanRank_CNV")
submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet Mean EXP Model", 
                       trainedModel= myMetabricENetCoxModelMeanRank_EXP,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelMeanRank_EXP.R", 
                       classFileName="metabricENetCoxModelMeanRank_EXP")
submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet IQM EXP CNV Model", 
                       trainedModel= myMetabricENetCoxModelIQMRank,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank.R", 
                       classFileName="metabricENetCoxModelIQMRank")
submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet IQM CNV Model", 
                       trainedModel= myMetabricENetCoxModelIQMRank_CNV,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank_CNV.R", 
                       classFileName="metabricENetCoxModelIQMRank_CNV")
submitCompetitionModel(modelName = "In Sock Parallel4 Rank M1 Enet IQM EXP Model", 
                       trainedModel= myMetabricENetCoxModelIQMRank_EXP,
                       rFiles="/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/parallelModel/round4/M/metabricENetCoxModelIQMRank_EXP.R", 
                       classFileName="metabricENetCoxModelIQMRank_EXP")
