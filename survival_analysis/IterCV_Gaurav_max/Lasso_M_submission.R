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

modelClassFile = ("~/COMPBIO/trunk/users/jang/R5/myEnetCoxIterMaxModel.R")
source(modelClassFile)

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/IterCV_Gaurav_max/ML/metabricENetCoxModelCancerCensusM.R")
modelClassFile2 = ("~/COMPBIO/trunk/users/jang/survival_analysis/IterCV_Gaurav_max/ML/metabricENetCoxModelMeanMarginalAssociationM.R")
modelClassFile3 = ("~/COMPBIO/trunk/users/jang/survival_analysis/IterCV_Gaurav_max/ML/metabricENetCoxModelMeanMetabricClusteringM.R")
modelClassFile4 = ("~/COMPBIO/trunk/users/jang/survival_analysis/IterCV_Gaurav_max/ML/metabricENetCoxModelMeanTopvaringHigginsM.R")
modelClassFile5 = ("~/COMPBIO/trunk/users/jang/survival_analysis/IterCV_Gaurav_max/ML/metabricENetCoxModelMeanTopvaringM.R")

source(modelClassFile1)
source(modelClassFile2)
source(modelClassFile3)
source(modelClassFile4)
source(modelClassFile5)

###################################################
### step 4: trainModel
###################################################
# Lasso Grid Setting
alphas = 1
lambdas = createENetTuneGrid(alphas = 1)[,2]
lambdas <- exp(seq(-5, 2, length = 50))

CancerCensusM <- metabricENetCoxModelCancerCensusM$new()
CancerCensusM$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions1 <- CancerCensusM$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationM <- metabricENetCoxModelMeanMarginalAssociationM$new()
MarginalAssociationM$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,iterNum = 20)
trainPredictions2 <- MarginalAssociationM$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringM <- metabricENetCoxModelMeanMetabricClusteringM$new()
MetabricClusteringM$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions3 <- MetabricClusteringM$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsM <- metabricENetCoxModelMeanTopvaringHigginsM$new()
TopvaringHigginsM$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions4 <- TopvaringHigginsM$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringM <- metabricENetCoxModelMeanTopvaringM$new()
TopvaringM$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions5 <- TopvaringM$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)


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
myModelName1 = "InSock Lasso max M CancerCensus" 
myModelName2 = "InSock Lasso max M MarginalAssociation" 
myModelName3 = "InSock Lasso max M MetabricClustering" 
myModelName4 = "InSock Lasso max M TopVaringHiggins" 
myModelName5 = "InSock Lasso max M TopVaring" 

submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusM,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationM,rFiles=list(modelClassFile2,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringM,rFiles=list(modelClassFile3,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsM,rFiles=list(modelClassFile4,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringM,rFiles=list(modelClassFile5,modelClassFile), parentDatasetId = "syn308537")

