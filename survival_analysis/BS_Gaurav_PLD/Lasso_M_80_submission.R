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

bsClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/BS_Gaurav_PLD/bsCoefM.R")

source(bsClassFile1)

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/BS_Gaurav_PLD/metabricENetCoxModelMprior.R")
source(modelClassFile1)

###################################################
### step 4: trainModel
###################################################
# Lasso Grid Setting

alphas = 1
lambdas <- exp(seq(-5, 2, length = 100))

betaFind1<-bsCoefM$new()
betaFind1$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "cancer_census",numBootstrap = 100)

betaFind2<-bsCoefM$new()
betaFind2$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "marginal",numBootstrap = 100)

betaFind3<-bsCoefM$new()
betaFind3$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "metabric",numBootstrap = 100)

betaFind4<-bsCoefM$new()
betaFind4$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvaringhiggins",numBootstrap = 100)

betaFind5<-bsCoefM$new()
betaFind5$customCoef(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvarings",numBootstrap = 100)

beta<- names(betaFind$coefBS)[which((betaFind$coefBS>=70)==TRUE)]

CancerCensusMprior <- metabricENetCoxModelMprior$new()
CancerCensusMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "cancer_census",numBootstrap = 100)
#trainPredictions1 <- CancerCensusMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MarginalAssociationMprior <- metabricENetCoxModelMprior$new()
MarginalAssociationMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "marginal",numBootstrap = 100)
# trainPredictions2 <- MarginalAssociationMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

MetabricClusteringMprior <- metabricENetCoxModelMprior$new()
MetabricClusteringMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "metabric",numBootstrap = 100)
# trainPredictions3 <- MetabricClusteringMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringHigginsMprior <- metabricENetCoxModelMprior$new()
TopvaringHigginsMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvaringhiggins",numBootstrap = 100)
# trainPredictions4 <- TopvaringHigginsMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

TopvaringMprior <- metabricENetCoxModelMprior$new()
TopvaringMprior$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas,priorName = "topvarings",numBootstrap = 100)
# trainPredictions5 <- TopvaringMprior$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

# 
# ###################################################
# ### step 5: computeTrainCIndex
# ###################################################
# trainPerformance1 <- SurvivalModelPerformance$new(as.numeric(trainPredictions1), trainingData$clinicalSurvData)
# trainPerformance2 <- SurvivalModelPerformance$new(as.numeric(trainPredictions2), trainingData$clinicalSurvData)
# trainPerformance3 <- SurvivalModelPerformance$new(as.numeric(trainPredictions3), trainingData$clinicalSurvData)
# trainPerformance4 <- SurvivalModelPerformance$new(as.numeric(trainPredictions4), trainingData$clinicalSurvData)
# trainPerformance5 <- SurvivalModelPerformance$new(as.numeric(trainPredictions5), trainingData$clinicalSurvData)
# 
# print(trainPerformance1$getExactConcordanceIndex())
# print(trainPerformance2$getExactConcordanceIndex())
# print(trainPerformance3$getExactConcordanceIndex())
# print(trainPerformance4$getExactConcordanceIndex())
# print(trainPerformance5$getExactConcordanceIndex())
# 
# 
# 
# ###################################################
# ### step 6: submitModel
# ###################################################
# myModelName1 = "InSock Lasso BS80 PLD M CancerCensus" 
# myModelName2 = "InSock Lasso BS80 PLD M MarginalAssociation" 
# myModelName3 = "InSock Lasso BS80 PLD M MetabricClustering" 
# myModelName4 = "InSock Lasso BS80 PLD M TopVaringHiggins" 
# myModelName5 = "InSock Lasso BS80 PLD M TopVaring" 
# 
# submitCompetitionModel(modelName = myModelName1, trainedModel=CancerCensusMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
# submitCompetitionModel(modelName = myModelName2, trainedModel=MarginalAssociationMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
# submitCompetitionModel(modelName = myModelName3, trainedModel=MetabricClusteringMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
# submitCompetitionModel(modelName = myModelName4, trainedModel=TopvaringHigginsMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
# submitCompetitionModel(modelName = myModelName5, trainedModel=TopvaringMprior,rFiles=list(modelClassFile1,modelClassFile), parentDatasetId = "syn308537")
# 
