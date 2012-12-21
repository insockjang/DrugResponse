### R code for predictive Modeling with Bayesian Lasso Regression comparing with Lasso
## RMSE is used to compare predictive performance
## Running Bayesian Lasso regression using {monomvn} package  
###################################################
### step 01: install package and load library
###################################################
# from CRAN
Install.packages("monomvn")

library(predictiveModeling)
library(affy)
library(monomvn)

# for preprocessing data (pruning missing value sample or uninformative features)
source("missingRow.R")
source("missingCol.R")

# Change the values of these variables
# starting with all data aggregated (cnv, mut, and expr)
myName <- "blasso_demo"
myWorkingDirectory <- "."
setwd(myWorkingDirectory)


###################################################
### step 02: Login Synapse for loading data and generate my own project
###################################################
library(synapseClient)

synapseLogin()

project <- Project(list(
  name=paste("Machine Learning Result (Bayesian Lasso Regression) - ", myName)
  ))
project <- createEntity(project)

analysis <- Analysis(list(
  name="blasso with All CCLE data",
  description="Bayesian Lasso Regression method run upon CCLE Data internal use only(not public)",
  parentId=propertyValue(project, "id")
	))
analysis <- createEntity(analysis)

dataset <- Dataset(list(
  name="Analysis Plots",
  parentId=propertyValue(project, "id")
  ))
dataset <- createEntity(dataset)

###################################################
#### Load Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "48339"     # rbin: "6067"
copyLayer <- loadEntity(id_copyLayer)
ds_copy_ccle <-exprs(copyLayer$objects$copySet)


id_oncomapLayer <- "48341"  # rbin "6069"
oncomapLayer <- loadEntity(id_oncomapLayer)
ds_oncomap_ccle <- exprs(oncomapLayer$objects$oncomapSet)

id_exprLayer <- "48344"     # rbin "6065"
exprLayer <- loadEntity(id_exprLayer)
ds_expr_ccle <- exprs(exprLayer$objects$exprSet)

###################################################
### Load Response Data
###################################################

id_drugLayer <- "48359" # sanger drug id "6084"
drugLayer <- loadEntity(id_drugLayer)
ds_drug_ccle <- exprs(drugLayer$objects$responseSet)


###################################################
### preprocessing data log transform
###################################################
# log transformation for sanger drug response data 
# need to compare log2 or log10 transformation vs. data preprocessing with add 1 
# meaningless values (~0 ) work as outlier or affect too much when training

ds_drug_ccle <- log10(ds_drug_ccle)

# need to be tested 
# ds_drug_ccle_1_log10 <- log10(ds_drug_ccle +1)
# ds_drug_ccle_log2 <- log2(ds_drug_ccle)
# ds_drug_ccle_1_log2 <- log2(ds_drug_ccle +1)


###################################################
### Aggregate all datasets (CNV + Mutation + Expression)
###################################################

ds_features_ccle <- createAggregateFeatureDataSet(list(copy = ds_copy_ccle, 
                                                        mut = ds_oncomap_ccle,
                                                        expr = ds_expr_ccle))
checkEquals(nrow(ds_features_ccle), nrow(ds_copy_ccle) + nrow(ds_oncomap_ccle) + nrow(ds_expr_ccle))

###################################################
### Data sample matching and pruning process for missing value samples or uninformative samples
###################################################
dataSets_ccle_Features_Response <- createFeatureAndResponseDataList(ds_features_ccle, 
                                                                           ds_drug_ccle)
ls(dataSets_ccle_Features_Response)
checkEquals(colnames(dataSets_ccle_Features_Response$featureData), 
            colnames(dataSets_ccle_Features_Response$responseData))

## To prune some features from feature space because of missing value problem 
mRow <- missingRow(dataSets_ccle_Features_Response$featureData) 

intResponse <- dataSets_ccle_Features_Response$responseData["PLX4720",,drop=FALSE]
mCol <- missingCol(intResponse)

featureData <- dataSets_ccle_Features_Response$featureData[-mRow,-mCol]
responseData <- intResponse[-mCol]

###################################################
## split the dataset for training and testing
## 2/3 portion will be chosen as traning set and the other will be assigned as testing set
###################################################

z <- sort(sample(1:length(responseData),floor(length(responseData) *2/3)))

mRow1 <- missingRow(featureData[,z]) 
mRow2 <- missingRow(featureData[,-z])

mRow <-union(mRow1,mRow2)

FeatureData_training <- featureData[-mRow,z]
ResponseData_training <- responseData[z]

FeatureData_testing <- featureData[-mRow,-z]
ResponseData_testing <- responseData[-z]

###################################################
### predictive Modeling core functions : bayesian lasso regression
###################################################

model_blasso1 <- blasso(t(FeatureData_training),ResponseData_training, T =1000, verb = 0)

# S <- summary(model_blasso)

BURNIN <- 200 # discard first 10% of MCMC sampling as a burning-in process

beta <- apply(model_blasso$beta[(BURNIN +1):nrow(model_blasso$beta),],2,mean)
beta <- as.numeric(beta)
intercept <- mean(model_blasso$mu[(BURNIN+1):length(model_blasso$mu)])
intercept <- as.numeric(intercept)

y_hat <-cbind(rep(1,dim(t(FeatureData_testing))[[1]]),t(FeatureData_testing)) %*% c(intercept,beta)

RMSE.blasso <-RMSE(y_hat,ResponseData_testing)


model_eNet <- fitPredictiveModel(t(FeatureData_training), 
                                 ResponseData_training, 
                                 method="glmnet", 
                                 tuneGrid=createENetTuneGrid(alphas=1)) # lasso
 
intercept_eNet <-model_eNet$finalModel$a0[length(model_eNet$finalModel$a0)]
y_hat2 <-cbind(rep(1,dim(t(FeatureData_testing[rownames(model_eNet$finalModel$beta),]))[[1]]),t(FeatureData_testing[rownames(model_eNet$finalModel$beta),])) %*% c(intercept_eNet,model_eNet$finalModel$beta[,ncol(model_eNet$finalModel$beta)])
RMSE.lasso <-RMSE(y_hat2,ResponseData_testing)

### NEED TO MODIFY From here to End ###############
### generate figure file
###################################################
coefs_blasso <- apply(model_blasso$beta,2,mean)
outputFileBayesianLasso <- 'PLX4720_BayesianLassoModel.jpg'
plotPredictiveModelHeatmap(coefs_blasso, FeatureData_training, ResponseData_training)

## predictive errors 
plot(ResponseData_testing-y_hat, col=2, pch=20)
points(ResponseData_testing-y_hat2, col=3, pch=18)

legend("topleft", c("blasso", "lasso"),col=c(2,3), pch=c(20,18))

###################################################
### store Result : need to implement efficiently
###################################################
blassoLayer <- Layer(list(
                            name="Bayesian Lasso Regression Results for PLX4720",
                            type="M", 
                            parentId=propertyValue(dataset, "id")))
blassoLayer <- addFile(blassoLayer, outputFileBayesianLasso)
storeEntity(blassoLayer)

step1 <- stopStep()
onWeb(step1)

propertyValue(step1, 'name') <- "Single run using bayesian lasso"
propertyValue(step1, 'description') <- "I found that ... looked very interesting due to ..."
step1 <- updateEntity(step1)
 
## bayesian lasso is computationally heavy, so let us discard cross validation part
#step2 <- startStep(analysis)
#propertyValue(step2, 'name') <- "Cross validation using elastic net"
#propertyValue(step2, 'input') <- propertyValue(step1, "input")
#step2 <- updateEntity(step2)
 
onWeb(analysis)

###################################################
### how to implement blasso into R class?
###################################################

setClass("BayesianLasso", 
    representation=representation(blassoModel="blasso"),
    contains="PredictiveModel",
    prototype = prototype(
			blassoModel = {
				nulllm <- list()
				class(nulllm) <- "blasso"
				nulllm
			}
	)
)

setMethod(
		f = "customTrain",
		signature = signature("BayesianLasso", "matrix", "numeric"),
		definition = function(method, featureData, responseData){
			message("in train_BayesianLasso")
			blassoModel <- blasso(responseData ~ featureData)
			names(blassoModel$coefficients)[2:length(blassoModel$coefficients)] <- colnames(featureData)
			method@blassoModel <- blassoModel
			method
		}
)

setMethod(
		f = "customPredict",
		signature = signature("BayesianLasso", "matrix"),
		definition = function(method, featureData){
			print("in predict_BayesianLasso")
			testCoefs <- method@blassoModel$coefficients
			testCoefs <- testCoefs[-1]
			testCoefs <- testCoefs[!is.na(testCoefs)]
			testCoefs <- as.matrix(testCoefs)
			testFeatures <- featureData[, row.names(testCoefs)]
			prediction <- testFeatures %*% testCoefs
			prediction <- prediction + method@blassoModel$coefficients[1]
			return(prediction)
		}
)

BayesianLasso <- function(){
  newObj <- new("BayesianLasso")
}


###################################################
### sessionInfo
###################################################
toLatex(sessionInfo())
