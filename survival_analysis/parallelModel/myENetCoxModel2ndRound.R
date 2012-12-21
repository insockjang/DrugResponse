library(predictiveModeling)
library(synapseClient)
library(survival)
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
fea<-t(unique(t(FEA))) # CNV data have lots of replicates

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

load("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/bestCombination.Rdata")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")

# make Parameter Grids
alphas <- c(10^(-10:-1), seq(0.2, 1, by = 0.1))
lambdas <- exp(seq(-10, 2, length = 30))  

######################################################################################################################################################
#######  Model 1. penalized CC with unpenalized best combination via IQM CC Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
myENetCoxModel_CCpenalizedCCunpenalizedIQM <- myEnetCoxModel$new()
myENetCoxModel_CCpenalizedCCunpenalizedIQM$customTrain(RES, survObj,alpha = alphas, lambda = lambdas, nfolds =5, penalty.factor = !best_IQM_combination[1,])
cvMyENetCoxModel_CCpenalizedCCunpenalizedIQM <- crossValidatePredictiveModel1(RES,survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5, penalty.factor = !best_IQM_combination[1,])

######################################################################################################################################################
#######  Model 2. penalized CC with unpenalized best combination via Mean CC Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
myENetCoxModel_CCpenalizedCCunpenalizedMean <- myEnetCoxModel$new()
myENetCoxModel_CCpenalizedCCunpenalizedMean$customTrain(RES, survObj,alpha = alphas, lambda = lambdas, nfolds =5, penalty.factor = !best_mean_combination[1,])
cvMyENetCoxModel_CCpenalizedCCunpenalizedMean <- crossValidatePredictiveModel1(RES,survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5, penalty.factor = !best_mean_combination[1,])

######################################################################################################################################################
#######  Model 3. penalized MF only Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
myENetCoxModel_MFonly <- myEnetCoxModel$new()
myENetCoxModel_MFonly$customTrain(FEA, survObj,alpha = alphas, lambda = lambdas, nfolds =5 )

cvMyENetCoxModel_MFonly <- crossValidatePredictiveModel1(FEA,survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5 )
######################################################################################################################################################
#######  Model 4. penalized CC+MF Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
myENetCoxModel_CCMFpenalized <- myEnetCoxModel$new()
myENetCoxModel_CCMFpenalized$customTrain(cbind(FEA,RES), survObj,alpha = alphas, lambda = lambdas, nfolds =5 )

cvMyENetCoxModel_CCMFpenalized <- crossValidatePredictiveModel1(cbind(FEA,RES),survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5 )

######################################################################################################################################################
#######  Model 5. unpenalized CC (best combination with IQM) + penalized MF Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
RRR<-RES[,which(best_IQM_combination[1,] !=0)]

myENetCoxModel_CCIQMunpenalizedMFpenalized <- myEnetCoxModel$new()
myENetCoxModel_CCIQMunpenalizedMFpenalized$customTrain(cbind(FEA,RRR), survObj,alpha = alphas, lambda = lambdas, nfolds =5,penalty.factor = c(rep(1,ncol(FEA)),rep(0,ncol(RRR))))
cvMyENetCoxModel_CCIQMunpenalizedMFpenalized <- crossValidatePredictiveModel1(cbind(FEA,RRR),survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5,penalty.factor=c(rep(1,ncol(FEA)),rep(0,ncol(RRR))))

######################################################################################################################################################
#######  Model 6. unpenalized CC (best combination with mean) + penalized MF Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
RRR<-RES[,which(best_mean_combination[1,] !=0)]

myENetCoxModel_CCMeanunpenalizedMFpenalized <- myEnetCoxModel$new()
myENetCoxModel_CCMeanunpenalizedMFpenalized$customTrain(cbind(FEA,RRR), survObj,alpha = alphas, lambda = lambdas, nfolds =5,penalty.factor = c(rep(1,ncol(FEA)),rep(0,ncol(RRR))))
cvMyENetCoxModel_CCMeanunpenalizedMFpenalized <- crossValidatePredictiveModel1(cbind(FEA,RRR),survObj,model = myEnetCoxModel$new(),numFolds =5, alpha = alphas, lambda= lambdas, nfolds =5,penalty.factor=c(rep(1,ncol(FEA)),rep(0,ncol(RRR))))

##############
##############
##############
## Analysis
##############
##############
##############

CV <- cvMyENetCoxModel_CCpenalizedCCunpenalizedIQM
cIndexTrain <-c()
cIndexTest <-c()

for(k in 1:5){
  cIndexTest<-c(cIndexTest,survConcordance(CV[[k]]$testObservations ~CV[[k]]$testPredictions)$concordance)
  cIndexTrain<-c(cIndexTrain,survConcordance(CV[[k]]$trainObservations ~CV[[k]]$trainPredictions)$concordance)
}
boxplot(cIndexTrain,cIndexTest)