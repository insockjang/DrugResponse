library(predictiveModeling)
library(synapseClient)
synapseLogin()

### need to source in this file since it is not in the predictiveModeling package
# source("/home/ijang/COMPBIO/trunk/users/jang/bootstrap/bootstrapPredictiveModel.R")

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

id_drugLayer <- "48337" # sanger drug response (ccle : "48359" but restricted access)
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$sangerADF

#data(demoData)

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
  
### reduce numBootstrap if you want to do faster testing
bootFeatureSelection <- bootstrapPredictiveModel(featureData_scaled,
                                                 responseData_scaled, 
                                                 numBootstrap= 100, 
                                                 model = CaretModel$new(modelType = "glmnet"), 
                                                 trControl = defaultTrainControl(), 
                                                 tuneGrid = createENetTuneGrid(alphas=1))

# Need function to select or give threshold to select significant features from result
plot(bootFeatureSelection)
significantFeatures <- names(bootFeatureSelection[which(bootFeatureSelection >= 9)]) # here I manually set the threshold 90% of bootNumbers
