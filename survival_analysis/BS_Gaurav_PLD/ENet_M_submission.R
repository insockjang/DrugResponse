### ENet Molecular Feature only 

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

modelClassFile = ("~/COMPBIO/trunk/users/jang/R5/myCoxModelBS.R")
source(modelClassFile)

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/BS_Gaurav_PLD/metabricENetCoxModelMprior.R")

source(modelClassFile1)

###################################################
### step 4: trainModel
###################################################
# ENet Grid Setting

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]
lambdas <- exp(seq(-5, 2, length = 100))

CancerCensusMprior <- metabricENetCoxModelMprior$new()
CancerCensusMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "cancer_census",numBootstrap = 100,threshold = 50)
trainPredictions1 <- CancerCensusMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationMprior <- metabricENetCoxModelMprior$new()
MarginalAssociationMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "marginal",numBootstrap = 100,threshold = 50)
trainPredictions2 <- MarginalAssociationMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringMprior <- metabricENetCoxModelMprior$new()
MetabricClusteringMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "metabric",numBootstrap = 100,threshold = 50)
trainPredictions3 <- MetabricClusteringMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsMprior <- metabricENetCoxModelMprior$new()
TopvaringHigginsMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvaringhiggins",numBootstrap = 100,threshold = 50)
trainPredictions4 <- TopvaringHigginsMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringMprior <- metabricENetCoxModelMprior$new()
TopvaringMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvarings",numBootstrap = 100,threshold = 50)
trainPredictions5 <- TopvaringMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)


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
myModelName1 = "InSock ENet BS PLD M CancerCensus" 
myModelName2 = "InSock ENet BS PLD M MarginalAssociation" 
myModelName3 = "InSock ENet BS PLD M MetabricClustering" 
myModelName4 = "InSock ENet BS PLD M TopVaringHiggins" 
myModelName5 = "InSock ENet BS PLD M TopVaring" 

submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")

