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
id_exprLayer <- "syn480759" 
layer_expr <- loadEntity(id_exprLayer)
data_expr <- layer_expr$objects$FoldChangeResponse

id_drugLayer <- "syn480806" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug<-layer_drug$objects$Phenotypes
data_drug<-as.matrix(adf_drug$BSEP.Inhibitor)
rownames(data_drug)<-adf_drug$SampleID

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])
                                        
filteredFeatureData  <- filteredData$featureData
filteredResponseData <- filteredData$responseData

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

resultsLasso<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=1, lambda = lambdas)
resultsENet<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=alphas, lambda = lambdas)
resultsRandomForest_100<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myRandomForestModel$new(), numFolds = 5,ntree = 100)
resultsRandomForest_1000<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myRandomForestModel$new(), numFolds = 5,ntree = 1000)
resultsRandomForest_10000<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myRandomForestModel$new(), numFolds = 5,ntree = 10000)

filename = paste("~/COMPBIO/trunk/users/jang/",Lara_BSEP,"/ROCR.Rdata",sep = "")
save(resultsLasso,resultsENet,resultsRandomForest,file = filename)

