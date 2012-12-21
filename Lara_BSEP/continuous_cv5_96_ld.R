### LARA BSEP classifier 
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
library(randomForest)
require(multicore)
require(doMC)
registerDoMC()

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical.R")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel1.R")

###################################################
#### Load BSEP Molecular Feature Data from Synapse
#### Categorical CV with 96 Hours
###################################################
data_expr <- train96ld
data_drug<-as.matrix(train96IC50)

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])

filteredFeatureData  <- filteredData$featureData
filteredResponseData <- filteredData$responseData

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

resultsLasso<-crossValidatePredictiveModel1(filteredFeatureData, filteredResponseData, model = myEnetModel$new(), alpha=1, lambda = lambdas)
resultsENet<-crossValidatePredictiveModel1(filteredFeatureData, filteredResponseData, model = myEnetModel$new(), alpha=alphas, lambda = lambdas)
resultsRandomForest_300<-crossValidatePredictiveModel1(filteredFeatureData, filteredResponseData, model = myRandomForestModel1$new(), numFolds = 5,ntree = 300)


catenateResults<-function(results){
  TestOb<-c()
  TestPr<-c()
  for(k in 1:length(results)){
    TestPr<-c(TestPr,results[[k]]$testPredictions)
    TestOb<-c(TestOb,results[[k]]$testObservations)  
  }
  return(list(cor = cor(TestPr,TestOb,use="complete"),
              predicted = TestPr,
              observed = TestOb))
}

resultsLasso<-catenateResults(resultsLasso)
resultsENet<-catenateResults(resultsENet)
resultsRandomForest_300<-catenateResults(resultsRandomForest_300)

filename = paste("~/COMPBIO/trunk/users/jang/Lara_BSEP/ROCR_continuous_cv5_96_ld.Rdata",sep = "")
save(resultsLasso,resultsENet,file = filename)
