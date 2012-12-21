library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

id_exprLayer <- "syn1532991" 
layer_expr <- loadEntity(id_exprLayer)
res2 <- layer_expr$objects$res2

id_headerLayer <- "syn1532987" 
layer_header <- loadEntity(id_headerLayer)
header <- layer_header$objects$header
# 

library(predictiveModeling)
library(ggplot2)
library(ROCR)
library(modeest)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/bootstrapPredictiveModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")

data_expr <- res2
data_drug<-as.numeric(as.matrix(header$DX))
names(data_drug)<-header$SAMPLE.ID

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,drop=FALSE])
filteredFeatureData  <- scale(filteredData$featureData)
filteredResponseData <- (filteredData$responseData)

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

resultsLasso<-bootstrapPredictiveModel(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), numBootstrap=100, alpha=1, lambda = lambdas)
resultsENet<-bootstrapPredictiveModel(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)

save(resultsLasso,file = "~/COMPBIO/trunk/users/jang/AD/Control_ADlike/bs100Lasso.Rdata")
save(resultsENet,file = "~/COMPBIO/trunk/users/jang/AD/Control_ADlike/bs100ENet.Rdata")
