rm(list = ls())
load("/home/cgaiteri/Amanda_A2Z/insock_brain.rdata")

library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

id_exprLayer <- "syn1532991" 
layer_expr <- loadEntity(id_exprLayer)
res2 <- layer_expr$objects$res2

id_headerLayer <- "syn1532987" 
layer_header <- loadEntity(id_headerLayer)
header <- layer_header$objects$header
# 
# ### check the match with colname and header order
# colnames(res2)[1:10]
# header$SAMPLE.ID[1:10]
# 
# ### split data into control vs. AD-like
# BinaryPhenotype<-header$DX
# cont <- which(BinaryPhenotype==0)
# 
# res2CONT<-res2[,cont]
# res2AD<-res2[,-cont]
# 
# ### check two datasets's relationship
# corcor<-c()
# for(k in 1:10){
#   A<-sample(nrow(res2),1000)
#   
#   corCONT<-cor(t(res2CONT[A,]))
#   corAD<-cor(t(res2AD[A,]))
#   
#   corcor<-c(corcor,cor(corCONT,corAD))
# }
# hist(corcor)
# # this looks like random
# 

### Do some classification works for 5 fold cross validation

library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical2.R")
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

resultsLasso<-crossValidatePredictiveModel_categorical2(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=1, lambda = lambdas)
resultsENet<-crossValidatePredictiveModel_categorical2(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=alphas, lambda = lambdas)

save(resultsLasso,file = "~/COMPBIO/trunk/users/jang/AD/Control_ADlike/cv5Lasso.Rdata")
save(resultsENet,file = "~/COMPBIO/trunk/users/jang/AD/Control_ADlike/cv5ENet.Rdata")
