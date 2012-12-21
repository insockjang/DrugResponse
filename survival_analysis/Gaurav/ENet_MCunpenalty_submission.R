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

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_unpenalty/metabricENetCoxModelCancerCensusMC_unpenalty.R")
modelClassFile2 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_unpenalty/metabricENetCoxModelMeanMarginalAssociationMC_unpenalty.R")
modelClassFile3 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_unpenalty/metabricENetCoxModelMeanMetabricClusteringMC_unpenalty.R")
modelClassFile4 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_unpenalty/metabricENetCoxModelMeanTopvaringHigginsMC_unpenalty.R")
modelClassFile5 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Gaurav/MC_unpenalty/metabricENetCoxModelMeanTopvaringMC_unpenalty.R")

source(modelClassFile1)
source(modelClassFile2)
source(modelClassFile3)
source(modelClassFile4)
source(modelClassFile5)

###################################################
### step 4: trainModel
###################################################
# ENet Grid Setting
alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

CancerCensusMC_unpenalty <- metabricENetCoxModelCancerCensusMC_unpenalty$new()
CancerCensusMC_unpenalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions1 <- CancerCensusMC_unpenalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationMC_unpenalty <- metabricENetCoxModelMeanMarginalAssociationMC_unpenalty$new()
MarginalAssociationMC_unpenalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions2 <- MarginalAssociationMC_unpenalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringMC_unpenalty <- metabricENetCoxModelMeanMetabricClusteringMC_unpenalty$new()
MetabricClusteringMC_unpenalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions3 <- MetabricClusteringMC_unpenalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsMC_unpenalty <- metabricENetCoxModelMeanTopvaringHigginsMC_unpenalty$new()
TopvaringHigginsMC_unpenalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions4 <- TopvaringHigginsMC_unpenalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringMC_unpenalty <- metabricENetCoxModelMeanTopvaringMC_unpenalty$new()
TopvaringMC_unpenalty$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions5 <- TopvaringMC_unpenalty$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)


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
myModelName1 = "InSock ENet MC_unpenalty CancerCensus" 
myModelName2 = "InSock ENet MC_unpenalty MarginalAssociation" 
myModelName3 = "InSock ENet MC_unpenalty MetabricClustering" 
myModelName4 = "InSock ENet MC_unpenalty TopVaringHiggins" 
myModelName5 = "InSock ENet MC_unpenalty TopVaring" 

submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusMC_unpenalty,rFiles=list(modelClassFile1), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationMC_unpenalty,rFiles=list(modelClassFile2), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringMC_unpenalty,rFiles=list(modelClassFile3), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsMC_unpenalty,rFiles=list(modelClassFile4), isPracticeModel=TRUE)
submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringMC_unpenalty,rFiles=list(modelClassFile5), isPracticeModel=TRUE)
