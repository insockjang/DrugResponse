### Lasso Molecular Feature only 

###################################################
### step 1: loadLibraries
###################################################
library(predictiveModeling)
library(BCC)
library(survival)
library(survcomp)
library(MASS)
library(rms)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")


###################################################ac
### step 2: loadData
###################################################
# synapseLogin() ### not required if configured for automatic login
trainingData <- loadMetabricTrainingData()


###################################################
### step 3: call predefined Models' classFile
###################################################

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_penalty/metabricENetCoxModelCancerCensusMC_penalty.R")
modelClassFile2 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_penalty/metabricENetCoxModelMeanMarginalAssociationMC_penalty.R")
modelClassFile3 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_penalty/metabricENetCoxModelMeanMetabricClusteringMC_penalty.R")
modelClassFile4 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_penalty/metabricENetCoxModelMeanTopvaringHigginsMC_penalty.R")
modelClassFile5 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_penalty/metabricENetCoxModelMeanTopvaringMC_penalty.R")

source(modelClassFile1)
source(modelClassFile2)
source(modelClassFile3)
source(modelClassFile4)
source(modelClassFile5)



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

CancerCensusMC_penalty <- metabricENetCoxModelCancerCensusMC_penalty$new()
CancerCensusMC_penalty$customTrain(exprData,copyData,clinicalFeaturesData,clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions1 <- CancerCensusMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

###################################################
### step 4: trainModel
###################################################
# Lasso Grid Setting
alphas = 1
lambdas = createENetTuneGrid(alphas = 1)[,2]

CancerCensusMC_penalty <- metabricENetCoxModelCancerCensusMC_penalty$new()
CancerCensusMC_penalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions1 <- CancerCensusMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationMC_penalty <- metabricENetCoxModelMeanMarginalAssociationMC_penalty$new()
MarginalAssociationMC_penalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions2 <- MarginalAssociationMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringMC_penalty <- metabricENetCoxModelMeanMetabricClusteringMC_penalty$new()
MetabricClusteringMC_penalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions3 <- MetabricClusteringMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsMC_penalty <- metabricENetCoxModelMeanTopvaringHigginsMC_penalty$new()
TopvaringHigginsMC_penalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions4 <- TopvaringHigginsMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringMC_penalty <- metabricENetCoxModelMeanTopvaringMC_penalty$new()
TopvaringMC_penalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions5 <- TopvaringMC_penalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)


###################################################
### step 5: computeTrainCIndex
###################################################
trainPerformance1 <- SurvivalModelPerformance$new(as.numeric(trainPredictions1), trainingData$clinicalSurvData)
trainPerformance2 <- SurvivalModelPerformance$new(as.numeric(trainPredictions2), trainingData$clinicalSurvData)
trainPerformance3 <- SurvivalModelPerformance$new(as.numeric(trainPredictions3), trainingData$clinicalSurvData)
trainPerformance4 <- SurvivalModelPerformance$new(as.numeric(trainPredictions4), trainingData$clinicalSurvData)
trainPerformance5 <- SurvivalModelPerformance$new(as.numeric(trainPredictions5), trainingData$clinicalSurvData)

print(trainPerformance1$getExactConcordanceIndex())
print(trainPerformance2$getExactConcordanceIndex())
print(trainPerformance3$getExactConcordanceIndex())
print(trainPerformance4$getExactConcordanceIndex())
print(trainPerformance5$getExactConcordanceIndex())



###################################################
### step 6: submitModel
###################################################
myModelName1 = "InSock Lasso MC_penalty CancerCensus" 
myModelName2 = "InSock Lasso MC_penalty MarginalAssociation" 
myModelName3 = "InSock Lasso MC_penalty MetabricClustering" 
myModelName4 = "InSock Lasso MC_penalty TopVaringHiggins" 
myModelName5 = "InSock Lasso MC_penalty TopVaring" 

submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusMC_penalty,rFiles=list(modelClassFile1), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationMC_penalty,rFiles=list(modelClassFile2), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringMC_penalty,rFiles=list(modelClassFile3), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsMC_penalty,rFiles=list(modelClassFile4), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringMC_penalty,rFiles=list(modelClassFile5), isPracticeModel=TRUE)
