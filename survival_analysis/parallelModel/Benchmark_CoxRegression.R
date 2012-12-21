library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@") ### not required if configured for automatic login

idExpressionLayer <- "139167"
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idCopyLayer <- "139169"
copyLayer <- loadEntity(idCopyLayer)
copyData <- copyLayer$objects[[1]]

idClinicalLayer <- "139171"
clinicalLayer <- loadEntity(idClinicalLayer)
clinicalData <- clinicalLayer$objects[[1]]
clinicalData <- clinicalData@data
clinicalData <- clinicalData[!is.na(clinicalData$survDeath), ]
#### prepare feature data for predictive modeling by transposing the matrix to have samples on the rows and features on the columns and scaling the columns
featureData <- scale(t(createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))))

featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "columns")

dataSets_expr_clinical <- createFeatureAndResponseDataList(featureData_filtered, clinicalData)

FEA <- dataSets_expr_clinical$featureData
RES <- dataSets_expr_clinical$responseData


# clinical covariate numeric mapping
res<-matrix(0,ncol = 14,nrow = dim(RES)[1])
rownames(res)<-rownames(RES)
colnames(res)<-colnames(RES)[1:14]

res[,"Site"] <- as.numeric(RES[,"Site"])
res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
res[,"grade"] <- as.numeric(RES[,"grade"])
res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
res[,"chemo"] <- as.numeric(RES[,"chemo"])
res[,"hormone"] <- as.numeric(RES[,"hormone"])
res[,"radiation"] <- as.numeric(RES[,"radiation"])
res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
res[,"tripleNegative"] <- ifelse(RES[,"tripleNegative"]=="pos", 1,0)
res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                            ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                   ifelse(RES[,"histology"]== "Medullary",3,
                                          ifelse(RES[,"histology"]== "MixedHistology",4,5))))
                     
rm(RES)
RES<-res

survObj <- Surv(dataSets_expr_clinical$responseData[,"survYears"], dataSets_expr_clinical$responseData[,"survDeath"])

load("/Volumes/ijang/Norway_Clinical/CoxModel/bestSelectionCC/bestCombination.Rdata")

######################################################################################################################################################
#######  Benchmark. all CC only Model : Cox Regression #################################################################################################
######################################################################################################################################################
# The number of folds for cross validation
numFolds = 5
foldIndices <- createFolds(RES[,1],k=numFolds)
RRR<-data.frame(res)

concordanceIQM <- foreach (fold = foldIndices) %do%{
  ccTrain <-  RRR[-fold,]
  ccTest  <-  RRR[fold,]
  
  survObjTrain  <- survObj[-fold]
  survObjTest   <- survObj[fold]
  
  fit <- coxph(survObjTrain ~ ., ccTrain)
  
  # predicted hazard ratio
  coxTrain <- survConcordance(survObjTrain ~ predict(fit))$concordance
  coxTest  <- survConcordance(survObjTest ~ predict(fit, ccTest))$concordance
  
  return(list(cTrain=coxTrain,cTest=coxTest))
}

cIndexTrain<-foreach(k= 1:numFolds) %do% concordanceIQM[[k]]$cTrain
trainCIndex <- do.call("c",cIndexTrain)
cIndexTest<-foreach(k= 1:numFolds) %do% concordanceIQM[[k]]$cTest
testCIndex <- do.call("c",cIndexTest)
mean(testCIndex)
######################################################################################################################################################
#######  Model 1. the best combination via IQM CC only Model : Cox Regression #################################################################################################
######################################################################################################################################################
# The number of folds for cross validation

concordanceIQM <- foreach (fold = foldIndices) %do%{
  ccTrain <-  RRR[-fold,which(best_IQM_combination[1,]!=0)]
  ccTest  <-  RRR[fold,which(best_IQM_combination[1,]!=0)]
  
  survObjTrain  <- survObj[-fold]
  survObjTest   <- survObj[fold]
  
  fit <- coxph(survObjTrain ~ ., ccTrain)
  
  # predicted hazard ratio
  coxTrain <- survConcordance(survObjTrain ~ predict(fit))$concordance
  coxTest  <- survConcordance(survObjTest ~ predict(fit, ccTest))$concordance
  
  return(list(cTrain=coxTrain,cTest=coxTest))
}

cIndexTrain<-foreach(k= 1:numFolds) %do% concordanceIQM[[k]]$cTrain
trainCIndex <- do.call("c",cIndexTrain)
cIndexTest<-foreach(k= 1:numFolds) %do% concordanceIQM[[k]]$cTest
testCIndex <- do.call("c",cIndexTest)
mean(testCIndex)
######################################################################################################################################################

######################################################################################################################################################
#######  Model 2. the best combination via mean CC only Model : Cox Regression #################################################################################################
######################################################################################################################################################

concordanceMean <- foreach (fold = foldIndices) %do%{
  ccTrain <-  RRR[-fold,which(best_mean_combination[1,]!=0)]
  ccTest  <-  RRR[fold,which(best_mean_combination[1,]!=0)]
  
  survObjTrain  <- survObj[-fold]
  survObjTest   <- survObj[fold]
  
  fit <- coxph(survObjTrain ~ ., ccTrain)
  
  # predicted hazard ratio
  coxTrain <- survConcordance(survObjTrain ~ predict(fit))$concordance
  coxTest  <- survConcordance(survObjTest ~ predict(fit, ccTest))$concordance
  
  return(list(cTrain=coxTrain,cTest=coxTest))
}

cIndexTrain<-foreach(k= 1:numFolds) %do% concordanceMean[[k]]$cTrain
trainCIndex <- do.call("c",cIndexTrain)
cIndexTest<-foreach(k= 1:numFolds) %do% concordanceMean[[k]]$cTest
testCIndex <- do.call("c",cIndexTest)
mean(testCIndex)
######################################################################################################################################################

# In most cases, model2 almost always slightly outperformed benchmark and model1. 
