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

modelClassFile1 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Insock/MC_unpenalty/MC_unpenalty_Exp.R")
modelClassFile2 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Insock/MC_unpenalty/MC_unpenalty_CNV.R")
modelClassFile3 = ("~/COMPBIO/trunk/users/jang/survival_analysis/Insock/MC_unpenalty/MC_unpenalty_ExpCNV.R")

source(modelClassFile1)
source(modelClassFile2)
source(modelClassFile3)

###################################################
### step 4: trainModel
###################################################
# Lasso Grid Setting
alphas = 1
lambdas = createENetTuneGrid(alphas = 1)[,2]
lambdas <- exp(seq(-5, 2, length = 50))

myMC_unpenalty_Exp <- MC_unpenalty_Exp$new()
myMC_unpenalty_Exp$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions1 <- myMC_unpenalty_Exp$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

myMC_unpenalty_CNV <- MC_unpenalty_CNV$new()
myMC_unpenalty_CNV$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions2 <- myMC_unpenalty_CNV$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)

myMC_unpenalty_ExpCNV <- MC_unpenalty_ExpCNV$new()
myMC_unpenalty_ExpCNV$customTrain(trainingData$exprData,trainingData$copyData,trainingData$clinicalFeaturesData,trainingData$clinicalSurvData, alpha = alphas,lambda = lambdas)
trainPredictions3 <- myMC_unpenalty_ExpCNV$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)




###################################################
### step 5: computeTrainCIndex
###################################################
trainPerformance1 <- SurvivalModelPerformance$new(as.numeric(trainPredictions1), trainingData$clinicalSurvData)
trainPerformance2 <- SurvivalModelPerformance$new(as.numeric(trainPredictions2), trainingData$clinicalSurvData)
trainPerformance3 <- SurvivalModelPerformance$new(as.numeric(trainPredictions3), trainingData$clinicalSurvData)

print(trainPerformance1$getExactConcordanceIndex())
print(trainPerformance2$getExactConcordanceIndex())
print(trainPerformance3$getExactConcordanceIndex())



###################################################
### step 6: submitModel
###################################################
myModelName1 = "InSock Lasso MC_penalty EXP" 
myModelName2 = "InSock Lasso MC_penalty CNV" 
myModelName3 = "InSock Lasso MC_penalty EXPCNV" 

submitCompetitionModel(modelName = myModelName1, trainedModel=myMC_unpenalty_Exp,rFiles=list(modelClassFile1), parentDatasetId = "syn311276")
submitCompetitionModel(modelName = myModelName2, trainedModel=myMC_unpenalty_CNV,rFiles=list(modelClassFile2), parentDatasetId = "syn311276")
submitCompetitionModel(modelName = myModelName3, trainedModel=myMC_unpenalty_ExpCNV,rFiles=list(modelClassFile3), parentDatasetId = "syn311276")

