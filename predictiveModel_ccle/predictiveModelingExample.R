library(predictiveModeling)
library(monomvn)



data(demoData)

ds_features_cn_mut_ccle <- createAggregateFeatureDataSet(list(copy = copySet, mut = oncomapSet))
dataSets_ccleFeaturesCnMut_sangerChems <- createFeatureAndResponseDataList(ds_features_cn_mut_ccle, sangerADF)

# Let us use the raw data, not the scaled(normalized data)
# you can find the rationale in the R code : whereDataNormalize.R
featureData_raw <- dataSets_ccleFeaturesCnMut_sangerChems$featureData
responseData_raw <- dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE]

A<-filterPredictiveModelData(t(featureData_raw),responseData_raw)
featureData_filtered <- t(A$featureData)
responseData_filtered <- A$responseData

featureData_filtered_norm <- t(scale(t(featureData_filtered)))
responseData_filtered_norm <- t(scale(t(responseData_filtered)))

source("/home/ijang/COMPBIO/trunk/users/jang/R5/myBlassoModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myBridgeModel.R")


# BLasso
predictiveModel_myBlasso <- myBlassoModel$new()
predictiveModel_myBlasso$train(t(featureData_raw), t(responseData_raw),mcmcNum = 5, burnIn = 1)
coefs_myBlasso <- predictiveModel_myBlasso$coefficients[2:length(predictiveModel_myBlasso$coefficients)]

# BRidge
predictiveModel_myBridge <- myBridgeModel$new()
predictiveModel_myBridge$train(t(featureData), t(responseData),mcmcNum = 5, burnIn = 1)
coefs_myBridge <- predictiveModel_myBridge$coefficients[2:length(predictiveModel_myBridge$coefficients)]

#################################### Start:  Caret predictive models
# Elastic Net
predictiveModel_eNet <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet$train(t(featureData_raw), t(responseData_raw), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet <- predictiveModel_eNet$rawCaretModel()
coefs_eNet <- caretModel_eNet$finalModel$beta[, ncol(caretModel_eNet$finalModel$beta)]


# ICR : independent component regression
predictiveModel_ICR <- CaretModel$new(modelType = "icr")
predictiveModel_ICR$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.n.comp = c(1,2,3,4)))
caretModel_ICR <- predictiveModel_ICR$rawCaretModel()
coefs_ICR <- caretMdoel_ICR$finalModel$coefficients[2:length(caretModel_ICR$finalModel$coefficients)]

# PCR : principle component regression
predictiveModel_PCR <- CaretModel$new(modelType = "pcr")
predictiveModel_PCR$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.ncomp = c(1,2,3,4)))
caretModel_PCR <- predictiveModel_PCR$rawCaretModel()
coefs_PCR <- caretModel_PCR$finalModel$coefficients[2:length(caretModel_PCR$finalModel$coefficients)]


# PLS : partial least squares
predictiveModel_PLS <- CaretModel$new(modelType = "pls")
predictiveModel_PLS$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.ncomp = c(1,2,3,4)))
caretModel_PLS <- predictiveModel_PLS$rawCaretModel()
coefs_PLS <- caretModel_PLS$finalModel$coefficients


# SVM : support vector machine
predictiveModel_SVM <- CaretModel$new(modelType = "svmLinear")
predictiveModel_SVM$train(t(featureData_raw), t(responseData_raw))
caretModel_SVM <- predictiveModel_SVM$rawCaretModel()
#coefs_SVM <- caretModel_SVM$finalModel@  


# LS : least square(linear regression)
predictiveModel_LS <- CaretModel$new(modelType = "lm")
predictiveModel_LS$train(t(featureData_raw), t(responseData_raw))
caretModel_LS <- predictiveModel_LS$rawCaretModel()
coefs_LS <- caretModel_LS$finalModel$coefficients[2:length(caretModel_LS$finalModel$coefficients)]

## RLS : robust least square(robust linear regression) : singular is not implemented here
#predictiveModel_RLS <- CaretModel$new(modelType = "rlm")
#predictiveModel_RLS$train(t(featureData), t(responseData))
#caretModel_RLS <- predictiveModel_RLS$rawCaretModel()
#coefs_RLS <- caretModel_RLS$finalModel$coefficients[2:length(caretModel_RLS$finalModel$coefficients)]

################################### End:  Caret predictive models


################################### START:  My R5 Classes for predictive models
predictiveModel_myLs <- myLsModel$new()
predictiveModel_myLs$train(t(featureData_raw), t(responseData_raw))
coefs_myLs <- predictiveModel_myLs$coefficients[2:length(predictiveModel_myLs$coefficients)]

predictiveModel_myEnet <- myEnetModel$new()
predictiveModel_myEnet$train(t(featureData_raw), t(responseData_raw), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_myEnet <- predictiveModel_myEnet$rawCaretModel()
coefs_myEnet <- caretModel_myEnet$finalModel$beta[, ncol(caretModel_myEnet$finalModel$beta)]

# ICR : independent component regression
predictiveModel_myIcr <- myIcrModel$new()
predictiveModel_myIcr$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.n.comp = c(1,2,3,4)))
model_myIcr <- predictiveModel_myIcr$rawCaretModel()
coefs_myIcr <- model_myIcr$finalModel$coefficients[2:length(caretModel_myIcr$finalModel$coefficients)]

# PCR : principle component regression
predictiveModel_myPcr <- myPcrModel$new()
predictiveModel_myPcr$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.ncomp = c(1,2,3,4)))
model_myPcr <- predictiveModel_myPcr$rawCaretModel()
coefs_myPcr <- model_myPcr$finalModel$coefficients[2:length(model_myPcr$finalModel$coefficients)]


# PLS : partial least squares
predictiveModel_myPls <- myPlsModel$new()
predictiveModel_myPls$train(t(featureData_raw), t(responseData_raw), tuneGrid=expand.grid(.ncomp = c(1,2,3,4)))
model_myPls <- predictiveModel_myPls$rawCaretModel()
coefs_myPls <- model_myPls$finalModel$coefficients


# SVM : support vector machine
predictiveModel_mySvm <- mySvmModel$new()
predictiveModel_mySvm$train(t(featureData_raw), t(responseData_raw),tuneGrid=expand.grid(.C = c(0.001,seq(0.05,0.95,by = 0.05),0.999)))
model_mySvm <- predictiveModel_mySvm$rawCaretModel()

################################### END :  My R5 Classes for predictive models

################################### START: cross-validated predictive models  using comparison between Caret vs. my customized R5 classes
# Elastic Net
cvResults_eNet <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=CaretModel$new(modelType="glmnet"), tuneGrid=createENetTuneGrid(alphas=1), numFolds=3)
cvResults_eNet$plotPredAndObs()

cvResults_myEnet <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=myEnetModel$new(), tuneGrid=createENetTuneGrid(alphas=1), numFolds=3)
cvResults_myEnet$plotPredAndObs()

# Least Square
cvResults_Ls <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=CaretModel$new(modelType="lm"), numFolds=3)
cvResults_Ls$plotPredAndObs()

cvResults_myLs <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=myLsModel$new(), numFolds=3)
cvResults_myLs$plotPredAndObs()

# Partial Least Square
cvResults_Pls <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=CaretModel$new(modelType="pls"),tuneGrid=expand.grid(.ncomp = c(1,2,3,4)), numFolds=3)
cvResults_Pls$plotPredAndObs()

cvResults_myPls <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=myPlsModel$new(), tuneGrid=expand.grid(.ncomp = c(1,2,3,4)), numFolds=3)
cvResults_myPls$plotPredAndObs()

# Support Vector Machine Regression
cvResults_Svm <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=CaretModel$new(modelType="svmLinear"), tuneGrid=expand.grid(.C = c(0.001,seq(0.05,0.95,by = 0.05),0.999)), numFolds=3)
cvResults_Svm$plotPredAndObs()

cvResults_mySvm <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=mySvmModel$new(), tuneGrid=expand.grid(.C = c(0.001,seq(0.05,0.95,by = 0.05),0.999)), numFolds=3)
cvResults_mySvm$plotPredAndObs()



## Warning : mcmcNum and burnIn should be specified before running !!!
# Bayesian Lasso
cvResults_myBlasso <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=myBlassoModel$new(), mcmcNum =5000, burnIn = 1000, numFolds=3)
cvResults_myBlasso$plotPredAndObs()                                                
                                                
# Bayesian Ridge
cvResults_myBridge <- crossValidatePredictiveModel(t(featureData_raw), t(responseData_raw), model=myBridgeModel$new(), mcmcNum =5000, burnIn = 1000, numFolds=3)
cvResults_myBridge$plotPredAndObs()

                                                
### Bootstarpping test
# Elastic Net
source("/home/ijang/COMPBIO/trunk/users/jang/bootstrap/bootstrapPredictiveModel.R")
bbbsResults_eNet <- bootstrapPredictiveModel(t(featureData_raw), t(responseData_raw), model=CaretModel$new(modelType="glmnet"), tuneGrid=createENetTuneGrid(alphas=1))
plot(bsResults_eNet)