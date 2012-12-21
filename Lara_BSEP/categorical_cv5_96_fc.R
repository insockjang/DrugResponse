### LARA BSEP classifier 
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical.R")
source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")

###################################################
#### Load BSEP Molecular Feature Data from Synapse ####
###################################################

data_expr <- train96fc
Name<-names(train96catBSEP)
data_drug<-as.numeric(as.matrix(train96catBSEP))
names(data_drug)<-Name
# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])

filteredFeatureData  <- filteredData$featureData
filteredResponseData <- filteredData$responseData


alphas = seq(0.1,1,by = 0.1)
lambdas = createENetTuneGrid(alphas = 1)[,2]

resultsLasso<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=1, lambda = lambdas)
resultsENet<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=alphas, lambda = lambdas)
resultsRandomForest_300<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myRandomForestModel$new(), numFolds = 5,ntree = 300)


filename = paste("~/COMPBIO/trunk/users/jang/Lara_BSEP/ROCR_categorical_cv5_96_fc.Rdata",sep = "")
save(resultsLasso,resultsENet,resultsRandomForest_300,file = filename)

