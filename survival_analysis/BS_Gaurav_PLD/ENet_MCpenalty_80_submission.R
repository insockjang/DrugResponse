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

modelClassFile = ("~/COMPBIO/trunk/users/jang/R5/myCoxModelBS.R")
source(modelClassFile)

bsClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/BS_Gaurav_PLD/bsCoefMCpenalty.R")

source(bsClassFile1)

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/BS_Gaurav_PLD/metabricENetCoxModelMC_penaltyprior.R")
source(modelClassFile1)

###################################################
### step 4: trainModel
###################################################
# Lasso Grid Setting

alphas =unique(createENetTuneGrid()[,1])
lambdas <- exp(seq(-5, 2, length = 100))

betaFind1<-bsCoefMCpenalty.R$new()
betaFind1$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "cancer_census",numBootstrap = 100)

betaFind2<-bsCoefMCpenalty.R$new()
betaFind2$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "marginal",numBootstrap = 100)

betaFind3<-bsCoefMCpenalty.R$new()
betaFind3$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "metabric",numBootstrap = 100)

betaFind4<-bsCoefMCpenalty.R$new()
betaFind4$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvaringhiggins",numBootstrap = 100)

betaFind5<-bsCoefMCpenalty.R$new()
betaFind5$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvarings",numBootstrap = 100)


CancerCensusMC_penaltyprior <- metabricENetCoxModelMC_penaltyprior$new()
CancerCensusMC_penaltyprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "cancer_census",numBootstrap = 100,threshold = 80)
trainPredictions1 <- CancerCensusMC_penaltyprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationMC_penaltyprior <- metabricENetCoxModelMC_penaltyprior$new()
MarginalAssociationMC_penaltyprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "marginal",numBootstrap = 100,threshold = 80)
trainPredictions2 <- MarginalAssociationMC_penaltyprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringMC_penaltyprior <- metabricENetCoxModelMC_penaltyprior$new()
MetabricClusteringMC_penaltyprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "metabric",numBootstrap = 100,threshold = 80)
trainPredictions3 <- MetabricClusteringMC_penaltyprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsMC_penaltyprior <- metabricENetCoxModelMC_penaltyprior$new()
TopvaringHigginsMC_penaltyprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvaringhiggins",numBootstrap = 100,threshold = 80)
trainPredictions4 <- TopvaringHigginsMC_penaltyprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringMC_penaltyprior <- metabricENetCoxModelMC_penaltyprior$new()
TopvaringMC_penaltyprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvarings",numBootstrap = 100,threshold = 80)
trainPredictions5 <- TopvaringMC_penaltyprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)


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
myModelName1 = "InSock ENet BS80 PLD MC_penalty CancerCensus" 
myModelName2 = "InSock ENet BS80 PLD MC_penalty MarginalAssociation" 
myModelName3 = "InSock ENet BS80 PLD MC_penalty MetabricClustering" 
myModelName4 = "InSock ENet BS80 PLD MC_penalty TopVaringHiggins" 
myModelName5 = "InSock ENet BS80 PLD MC_penalty TopVaring" 

submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusMC_penaltyprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationMC_penaltyprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringMC_penaltyprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsMC_penaltyprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringMC_penaltyprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
