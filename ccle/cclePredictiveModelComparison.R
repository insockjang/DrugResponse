library(predictiveModeling)
library(synapseClient)

synapseLogin()

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

foreach(k = 1:dim(drugResponse)[2]) %dopar% {
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,k])
  featureData_processed <- filteredData$featureData
  responseData_processed <- filteredData$responseData
  
  ## scale this data
  featureData_scaled <- scale(featureData_processed)
  responseData_scaled <- scale(responseData_processed)
  
  predictiveModel_eNet <- myEnetModel$new()
  predictiveModel_eNet$customTrain(featureData_scaled,responseData_scaled)
  predictedTrain<-predictiveModel_eNet$customPredict(featureData_scaled)
  
  fit1<-predictiveModel_eNet$rawModel()
  plot(fit1)
  coefficient <- predictiveModel_eNet$getCoefficients()
  
  # cross validation for training testing set 
}
########  2/18/2012 need to implement furthermore
  
## fit a model
predictiveModel_eNet <- fitPredictiveModel(t(featureData_scaled), t(responseData_scaled), method="glmnet", tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet <- rawCaretModel(predictiveModel_eNet)
coefs_eNet <- caretModel_eNet$finalModel$beta[, ncol(caretModel_eNet$finalModel$beta)]
plotPredictiveModelHeatmap(coefs_eNet, featureData_scaled, responseData_scaled)

## fit a model with custom elastic net class
predictiveModel_my_eNet <- fitPredictiveModel(t(featureData_scaled), t(responseData_scaled), method=new("MyElasticNetModel"))
coefs_my_eNet <- predictiveModel_my_eNet@model$model$finalModel$beta[, ncol(predictiveModel_my_eNet@model$model$finalModel$beta)]
plotPredictiveModelHeatmap(coefs_my_eNet, featureData_scaled, responseData_scaled)

## Independent Component Regression <-(faster, stable)<- Principal Component Regression
# cvResult_PCR <- crossValidatePredictiveModel(t(featureData_scaled),t(responseData_scaled), method=new("MyPCRModel"), numFolds=3)

cvResult_ICR <- crossValidatePredictiveModel(t(featureData_scaled),t(responseData_scaled), method=new("MyICRModel"), numFolds=3)
cvResult_PLS <- crossValidatePredictiveModel(t(featureData_scaled),t(responseData_scaled), method=new("MyPLSModel"), numFolds=3)
cvResults_eNet <- crossValidatePredictiveModel(t(featureData_scaled), t(responseData_scaled), method="glmnet", tuneGrid=createENetTuneGrid(alphas=1), numFolds=3)
cvResults_myENet <- crossValidatePredictiveModel(t(featureData_scaled), t(responseData_scaled), method=new("MyElasticNetModel"), numFolds = 3)

### test errors should be the same between elastic net and custom class implementation of elastic net
plot(getTestError(cvResults_eNet), getTestError(cvResults_myENet))

### and confirm this numerically
getR2(cvResults_eNet)
getR2(cvResults_myENet)
getR2(cvResults_PCR)
getR2(cvResults_ICR)
getR2(cvResults_PLR)