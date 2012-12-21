# This code is used for checking data normailzation
# in some R packages such as elasticNet (glmnet), lars or (b)lasso/(b)ridge, the normalization procedure is automatically and internally applied when given input data are not normalized
# This normalization procedure is only used for feature selection. 
# The trained coefficients might be from row data !!! WARNING !!! need to check when applying predictive models)
# Without being careful, your predictive performance comparison means NOTHING...

library(predictiveModeling)

data(demoData)

# This is for aggregate the feature datasets from copy number and mutation 
ds_features_cn_mut_ccle <- createAggregateFeatureDataSet(list(copy = copySet, mut = oncomapSet))
dataSets_ccleFeaturesCnMut_sangerChems <- createFeatureAndResponseDataList(ds_features_cn_mut_ccle, sangerADF)

# noprocessedData 00
rawFeatureData <- dataSets_ccleFeaturesCnMut_sangerChems$featureData
rawResponseData <- dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE]


# preprocessedData 01 : filtering out 'NA' rows and columns
processedData <- filterPredictiveModelData(t(dataSets_ccleFeaturesCnMut_sangerChems$featureData),t(dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE]))
featureData <- t(processedData$featureData)
responseData <- t(processedData$responseData)

# preprocessedData02 : prenormalized data before applying(inputing) them for predictive models  
featureData_scaled <- t(scale(t(dataSets_ccleFeaturesCnMut_sangerChems$featureData)))
responseData_scaled <- t(scale(t(dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE])))

# preprocessingData03 : two combining steps
# 1. filtering out 'NA' rows and columns
# 2. normalize the filtered-out data
featureData_scaled1 <- t(scale(t(featureData)))
responseData_scaled1 <- t(scale(t(responseData)))

# 1. We need to check if the selected features are the same when elastic net is applied with preprocessedData01 and preprocessedData03
# 2. In this case, the estimated beta(coefficients) should be different becasue the response and training dataset are different  
# 3. Be careful -> When you compare the caretModels, the parameters should be the same (alpha, and lambda)

#################################### Start:  Caret predictive models
# Elastic Net
predictiveModel_eNet <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet$train(t(rawFeatureData), t(rawResponseData), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet <- predictiveModel_eNet$rawCaretModel()
coefs_eNet <- caretModel_eNet$finalModel$beta[, ncol(caretModel_eNet$finalModel$beta)]

predictiveModel_eNet1 <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet1$train(t(featureData), t(responseData), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet1 <- predictiveModel_eNet1$rawCaretModel()
coefs_eNet1 <- caretModel_eNet1$finalModel$beta[, ncol(caretModel_eNet1$finalModel$beta)]

predictiveModel_eNet2 <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet2$train(t(featureData_scaled), t(responseData_scaled), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet2 <- predictiveModel_eNet2$rawCaretModel()
coefs_eNet2 <- caretModel_eNet2$finalModel$beta[, ncol(caretModel_eNet2$finalModel$beta)]

predictiveModel_eNet3 <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet3$train(t(featureData_scaled1), t(responseData_scaled1), tuneGrid=createENetTuneGrid(alphas=1))
caretModel_eNet3 <- predictiveModel_eNet3$rawCaretModel()
coefs_eNet3 <- caretModel_eNet3$finalModel$beta[, ncol(caretModel_eNet3$finalModel$beta)]

# The result from coefs_eNet(from rawData) and coefs_eNet1(from filteredData) shows exactly the same results.
# Check if the selected features are different 
setdiff(names(coefs_eNet),names(coefs_eNet1))
setdiff(names(coefs_eNet),names(coefs_eNet2))
setdiff(names(coefs_eNet),names(coefs_eNet3))
setdiff(names(coefs_eNet1),names(coefs_eNet2))
setdiff(names(coefs_eNet1),names(coefs_eNet3))
setdiff(names(coefs_eNet2),names(coefs_eNet3))

par(mfrow = c(1, 3))
plot(coefs_eNet,coefs_eNet1)
plot(coefs_eNet,coefs_eNet2)
plot(coefs_eNet1,coefs_eNet2)




# For supporting idea that we should not normalize the input data. 
# if we preprocess(normalize) data and apply them for predictive models, the performance such as RMSE drastically changed (it is obvious !) 
# LS : least square(linear regression)
predictiveModel_LS <- CaretModel$new(modelType = "lm")
predictiveModel_LS$train(t(rawFeatureData), t(rawResponseData))
caretModel_LS <- predictiveModel_LS$rawCaretModel()
coefs_LS <- caretModel_LS$finalModel$coefficients[2:length(caretModel_LS$finalModel$coefficients)]

predictiveModel_LS1 <- CaretModel$new(modelType = "lm")
predictiveModel_LS1$train(t(featureData), t(responseData))
caretModel_LS1 <- predictiveModel_LS1$rawCaretModel()
coefs_LS1 <- caretModel_LS1$finalModel$coefficients[2:length(caretModel_LS$finalModel$coefficients)]

predictiveModel_LS2 <- CaretModel$new(modelType = "lm")
predictiveModel_LS2$train(t(featureData_scaled), t(responseData_scaled))
caretModel_LS2 <- predictiveModel_LS2$rawCaretModel()
coefs_LS2 <- caretModel_LS2$finalModel$coefficients[2:length(caretModel_LS2$finalModel$coefficients)]

predictiveModel_LS3 <- CaretModel$new(modelType = "lm")
predictiveModel_LS3$train(t(featureData_scaled1), t(responseData_scaled1))
caretModel_LS3 <- predictiveModel_LS2$rawCaretModel()
coefs_LS3 <- caretModel_LS3$finalModel$coefficients[2:length(caretModel_LS3$finalModel$coefficients)]

# Check if the selected features are different 
setdiff(names(coefs_LS[which(!is.na(coefs_LS))]),names(coefs_LS[which(!is.na(coefs_LS1))]))
setdiff(names(coefs_LS[which(!is.na(coefs_LS))]),names(coefs_LS[which(!is.na(coefs_LS2))]))
setdiff(names(coefs_LS[which(!is.na(coefs_LS1))]),names(coefs_LS[which(!is.na(coefs_LS2))]))
                                          
par(mfrow = c(1, 3))
plot(coefs_LS,coefs_LS1)
plot(coefs_LS,coefs_LS2)
plot(coefs_LS1,coefs_LS2)

