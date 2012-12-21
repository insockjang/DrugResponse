library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/ccle/crossValidatePredictiveModel_test.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/ENetModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/cvENetModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myCvEnetModel.R")
source('/home/ijang/COMPBIO/trunk/users/jang/R5/gridAlphaENetResults.R')
source('/home/ijang/COMPBIO/trunk/users/jang/ccle/cvGridAlphaMyENetModel.R')

source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/ccle/crossValidatePredictiveModel_test.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/ENetModel.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/cvENetModel.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/myCvEnetModel.R")
source('/Volumes/ijang/COMPBIO/trunk/users/jang/R5/gridAlphaENetResults.R')
source('/Volumes/ijang/COMPBIO/trunk/users/jang/ccle/cvGridAlphaMyENetModel.R')
library(glmnet)
###################################################
#### Load Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "48339"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$copySet

id_oncomapLayer <- "48341"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$oncomapSet

id_exprLayer <- "48344" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$exprSet

###################################################
### Load Response Data
###################################################

id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

drugResponse <- dataSets_ccle$responseData

# for all drug response prediction
# createENetTuneGrid need modified in the argument of "len"

alphas <- c(10^(-5:-1), seq(0.2, 1, by = 0.1))
c1 <-createENetTuneGrid(alphas = 1)
d1<-expand.grid(alpha = alphas,lambda = c1[,2])
d<-list()
for (k in 1:length(unique(d1[,1]))){
  d[[k]]<-expand.grid(alpha =unique(d1[,1])[k],unique(d1[,2]))
  
}

# Predictive step step without cross validation 
# small feature set are artificially selected in order to save computational time
gridEnetResults <- foreach(k = 1:ncol(dataSets_ccle$responseData)) %dopar%{  
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,k])
  featureData_processed <- filteredData$featureData
  responseData_processed <- filteredData$responseData
  
  ## scale this data
  featureData_scaled <- scale(featureData_processed)
  responseData_scaled <- scale(responseData_processed)
  fff<-featureData_scaled[,1:500]
  return(gridAlphaENetResults(fff, responseData_scaled, model = myEnetModel, alphas = alphas, lambda = d[[1]][,2],nfolds = 5))
}

# Predictive step step with cross validation
# small feature set are artificially selected in order to save computational time
cvGridEnetResults <- foreach(k = 1:ncol(dataSets_ccle$responseData)) %dopar%{
  
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,k])
  featureData_processed <- filteredData$featureData
  responseData_processed <- filteredData$responseData
  
  ## scale this data
  featureData_scaled <- scale(featureData_processed)
  responseData_scaled <- scale(responseData_processed)
  fff<-featureData_scaled[,1:500]
  return(cvGridAlphaMyENetModel(fff, responseData_scaled, model = myEnetModel, alphas = alphas, lambda = d[[1]][,2],nfolds = 5))
}

# Predictive step with cross validation
# Full feature set is used for real research
cvGridEnetResults <- foreach(k = 1:ncol(dataSets_ccle$responseData)) %dopar%{
  
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,k])
  featureData_processed <- filteredData$featureData
  responseData_processed <- filteredData$responseData
  
  ## scale this data
  featureData_scaled <- scale(featureData_processed)
  responseData_scaled <- scale(responseData_processed)
  
  return(cvGridAlphaMyENetModel(featureData_scaled, responseData_scaled, model = myEnetModel, alphas = alphas, lambda = d[[1]][,2],nfolds = 5))
}

