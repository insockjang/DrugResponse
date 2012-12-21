library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)

synapseLogin() ### not required if configured for automatic login

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

source("/Volumes/ijang/COMPBIO/trunk/users/jang/bootstrap/bootstrapGrid.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/crossValidationParameterOptimizer_penalty.R")
source("/Volumes/ijang/COMPBIO/trunk/users/jang/R5/ENetModel.R")

# make Parameter Grids
alpha = seq(0.05,1,0.05)
lambda = createENetTuneGrid(alphas=1)[,2]

######################################################################################################################################################
#######  Model 1. 1st find optimal alpha and lambda by Cross Validation  #############################################################################
#######           2nd bootstrap to find more significant features within selected feature            #################################################
######################################################################################################################################################
cModel2 <- crossValidationParameterOptimizer_penalty$new()
cModel2$findOpt(model=ENetModel, RES, survObj, numFolds = 10,  alpha = alpha, lambda = lambda, family = "cox")
cModel2$opt.metricGrid
cModel2$opt.param

Model2$opt.metricGrid
Model2$opt.param


######################################################################################################################################################
#######  Model 2. 1st bootstrap to each (alpha, lambda) in Grid space to find significant features with threshold    ############################
#######           2nd Run coxph model with selected features     #############################################################################
######################################################################################################################################################
A<-bootstrapFeatureGrid(model = ENetModel,RES,survObj,numBootstrap = 100,alpha =alpha,lambda= lambda,family = "cox")

B<-foreach(k = 1:length(alpha)) %do% {A[[k]] >= 70}

# Rerun with ENet algorithm to estimate beta(coefficients)
alphaModel <- foreach(k = 1:length(alpha)) %dopar% {
  lambdaModel <- foreach(j = 1:length(lambda)) %dopar%{
    predictiveModelCox<-ENetModel$new()    
    if (sum(B[[k]][,j])==0)
      return(predictiveModelCox)
    else{
      predictiveModelCox<-ENetModel$new()
      predictiveModelCox$customTrain(RES,
                                     survObj,
                                     alpha = alpha[k],
                                     lambda = lambda[j], 
                                     family = "cox",
                                     penalty.factor = as.numeric(xor(rep(1,ncol(RES)),as.numeric(B[[k]][,j]))))
      return(predictiveModelCox)
    }
  }
  return(lambdaModel)
}


# Run with Coxph with dataframe 
alphaFit <-foreach(k=1:length(alpha)) %dopar% {
  lambdaFit<-foreach(j=1:length(lambda)) %dopar% {
    if (sum(B[[k]][,j])==0)
      fit<-NA
    else
      fit<-coxph(survObj~.,res[,B[[k]][,j]])
    return(fit)
  }
  return(lambdaFit)
}

numFolds = 10
foldIndices <- createFolds(res[,1],k=numFolds)
MatTrain <- list()
MatTest <- list()

for(i in 1:length(foldIndices)){
  cMatTrain<-matrix(0,nrow = length(alpha),ncol = length(lambda))
  cMatTest<-matrix(0,nrow = length(alpha),ncol = length(lambda))
  for(k in 1:length(alpha)){
    for(j in 1:length(lambda)) {
      if (sum(B[[k]][,j])<=1){
        cMatTrain[k,j]<-0
        cMatTest[k,j]<-0
      }
      else{
        fit<-coxph(survObj[-foldIndices[[i]]]~.,res[-foldIndices[[i]],B[[k]][,j]])
        cIndexTrain<-concordance.index(predict(fit),survObj[-foldIndices[[i]],1],survObj[-foldIndices[[i]],2])
        cMatTrain[k,j]<- cIndexTrain$c.index
        cIndexTest<-concordance.index(predict(fit,res[foldIndices[[i]],B[[k]][,j]]),survObj[foldIndices[[i]],1],survObj[foldIndices[[i]],2])
        cMatTest[k,j]<- cIndexTest$c.index
      }
    }
  }
  MatTrain[[i]]<-cMatTrain
  MatTest[[i]]<-cMatTest
}

IQM <- function(vec,na.rm = TRUE){
  a<-quantile(vec,na.rm = TRUE)
  b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
  return(mean(vec[b],na.rm = TRUE))
}

MAT_train <- matrix(0,nrow = length(alpha),ncol = length(lambda))
MAT_test <- matrix(0,nrow = length(alpha),ncol = length(lambda))

for(k in 1:length(alpha)){
  for(j in 1:length(lambda)) {
    ccTrain<-c()
    ccTest<-c()
    for(i in 1:numFolds){
      ccTrain<-c(ccTrain,MatTrain[[i]][k,j])
      ccTest<-c(ccTest,MatTest[[i]][k,j])
    }
    MAT_train[k,j] = IQM(ccTrain)
    MAT_test[k,j] = IQM(ccTest)
  }
}


# for checking the coefficients for (alpha, lambda) in Grid space
#C1<-foreach(k=1:length(lambda))%do% as.matrix(alphaModel[[10]][[k]]$getCoefficients())
#CC1 <-do.call("cbind",C1)
