### DEMO DATA
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
library(standGL)
### need to source in this file since it is not in the predictiveModeling package
source("/home/ijang/COMPBIO/trunk/users/jang/prior_info/groupIndexForm.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myStandGLModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/StandGLModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/CrossValidationParameterOptimizer.R")

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

filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,"PLX4720",drop=FALSE])

# to check if this demo works, I make the feature dataset be very small (~ 500 features)
# featureData_processed <- filteredData$featureData[,1:500]

##### select small number of features for testing to make the demo run faster.
featureData_processed <- filteredData$featureData
responseData_processed <- filteredData$responseData

## scale this data
featureData_scaled <- scale(featureData_processed)
responseData_scaled <- scale(responseData_processed)



groupIndex <- groupIndexForm(featureData_scaled)
lambda <-createENetTuneGrid(alphas= 1)

# intercept included
GLmodel<-myStandGLModel$new()

# without intercept
GLmodel<-StandGLModel$new()
GLmodel$customTrain(featureData_scaled,responseData_scaled,index = groupIndex, lambda = lambda[,2], alpha = 1)



GLmodel$customTrain(featureData_scaled,responseData_scaled,index = index, lambda = lambda[,2], alpha = 1)

lambda <-createENetTuneGrid(alphas= 1)
index1<-sample(500,500,replace = TRUE)
GLmodel<-StandGLModel$new()
GLmodel$customTrain(featureData_scaled[,1:500],responseData_scaled,index = sort(index1), lambda = lambda[,2], alpha = 1)

# Yes, this is working !!!
cv1<-CrossValidationParameterOptimizer$new()
cv1$findOptimalParameters(model = StandGLModel,FFF,RRR,numFolds = 3,tuneGrid = createENetTuneGrid(alphas=1),index = sort(index1))


cv<-CrossValidationParameterOptimizer$new()
cv$findOptimalParameters(model = StandGLModel,featureData_scaled,responseData_scaled,numFolds = 10,tuneGrid = createENetTuneGrid(alphas=1),index = groupIndex)


