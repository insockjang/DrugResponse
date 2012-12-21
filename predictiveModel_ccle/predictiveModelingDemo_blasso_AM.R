### R code for predictive Modeling with Bayesian Lasso Regression comparing with Lasso
## RMSE is used to compare predictive performance
## Running Bayesian Lasso regression using {monomvn} package  
###################################################
### step 01: install package and load library
###################################################

library(predictiveModeling)
library(affy)
library(monomvn)

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

heatmap(t(scale(t(ds_drug_ccle))))
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


featureData <- dataSets_ccle_Features_Response$featureData
responseData <- dataSets_ccle_Features_Response$responseData["AZD6244",,drop=FALSE]

model_eNet <- fitPredictiveModel(t(featureData),responseData, method="glmnet", tuneGrid=createENetTuneGrid(alphas=1)) # lasso

model_eNet <- crossValidatePredictiveModel(t(featureData),responseData, method="glmnet", tuneGrid=createENetTuneGrid(alphas=1)) # lasso

par(mfrow = c(1,2))
plot(model_eNet$predAndObs_train)
plot(model_eNet$predAndObs_test)
###################################################
### predictive Modeling core functions : bayesian lasso regression
###################################################


###################################################
### how to implement blasso into R class?
###################################################
setClass("BayesianLassoRegression",
         
         representation=representation(blassoModel="blasso"),
         
         contains="PredictiveModel",
         
         prototype = prototype(
           
                      blassoModel = {

                                nullblasso <- list()

                                class(nullblasso) <- "blasso"

                                nullblasso

                        }

        )

)

 

setMethod(

                f = "customTrain",

                signature = signature("BayesianLassoRegression", "matrix", "numeric"),

                definition = function(method, featureData, responseData){

                        message("in train_BayesianLassoRegression")

                        blassoModel <- blasso(featureData,responseData,T=10,verb=0)
                        
                        BURNIN <- 5
                        # discard first 10% of MCMC sampling as a burning-in process

                        beta <- apply(blassoModel$beta[(BURNIN +1):nrow(blassoModel$beta),],2,mean)
                        
                        names(beta) <- colnames(featureData)
          
                        intercept <- mean(blassoModel$mu[(BURNIN+1):length(blassoModel$mu)])

                        blassoModel$coef <- beta
                        blassoModel$intercept <- intercept
                      
                        method@blassoModel <- blassoModel

                        method

                }

)

 

setMethod(

                f = "customPredict",

                signature = signature("BayesianLassoRegression", "matrix"),

                definition = function(method, featureData){

                        print("in predict_linearRegresion")

                        testCoefs <- method@blassoModel$coef

                        testCoefs <- as.matrix(testCoefs)

                        testFeatures <- featureData[, row.names(testCoefs)]

                        prediction <- testFeatures %*% testCoefs

                        prediction <- prediction + method@blassoModel$intercept

                        return(prediction)

                }

)

 

BayesianLassoRegression <- function(){

 newObj <- new("BayesianLassoRegression")

}

 

 

###################################################

### code chunk number 23: predictiveModelingDemo.Rnw:401-405

###################################################

blasso_cv <- crossValidatePredictiveModel(t(FeatureData),

                                           ResponseData,

                                           method=BayesianLassoRegression(),

                                           numFolds=10)

###################################################
### sessionInfo
###################################################
toLatex(sessionInfo())
