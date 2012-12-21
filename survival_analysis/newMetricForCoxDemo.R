library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/findOptimalCoxModel.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/findOptimalCoxParam_DIndex.R")
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/crossValidateFindOptimalCoxModel.R")
synapseLogin("in.sock.jang@sagebase.org")

# Load Data from Synapse
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

# Survival Object which is used for Response data in Cox Model
survObj <- Surv(dataSets_expr_clinical$responseData[,"survYears"], dataSets_expr_clinical$responseData[,"survDeath"])

# how many cross validate folds?
numFolds = 5
######################################################################################################################################################
#######  unpenalized CC only Model : Cox Regression #################################################################################################
######################################################################################################################################################
foldIndices <- createFolds(RES[,1],k=numFolds)

RRR<-data.frame(res)
coxTrain<-list()
coxTest<-list()
for (k in 1:5){
  resTrain<-RRR[-foldIndices[[k]],]
  resTest<-RRR[foldIndices[[k]],]
  
  survObjTrain <-Surv(dataSets_expr_clinical$responseData[-foldIndices[[k]],"survYears"], dataSets_expr_clinical$responseData[-foldIndices[[k]],"survDeath"])
  survObjTest <-Surv(dataSets_expr_clinical$responseData[foldIndices[[k]],"survYears"], dataSets_expr_clinical$responseData[foldIndices[[k]],"survDeath"])
  fit<-coxph(survObjTrain ~ Site + ageDiagnosis + lymphnodes + grade + histology + tumorSizeCM + chemo + hormone + radiation + HER2 +ER + PR + ERPR + tripleNegative, resTrain)
  #fit<-coxph(survObjTrain ~ ageDiagnosis + lymphnodes + grade + tumorSizeCM + HER2 + ERPR, resTrain)
  coxTrain[[k]] <- predict(fit)
  coxTest[[k]] <- predict(fit, resTest)
}
# by default, 5 fold cross validation
cindexTrain<-c()
cindexTest<-c()
dindexTrain<-c()
dindexTest<-c()
for (i in 1:numFolds){
  cTrain<-concordance.index(x=coxTrain[[i]], surv.time=survObj[-foldIndices[[i]],1], surv.event=survObj[-foldIndices[[i]],2], na.rm=TRUE, alpha= .05, method="noether")
  cindexTrain <- c(cindexTrain,cTrain$c.index)
  cTest<-concordance.index(x=coxTest[[i]], surv.time=survObj[foldIndices[[i]],1], surv.event=survObj[foldIndices[[i]],2], na.rm=TRUE, alpha= .05, method="noether")
  cindexTest <- c(cindexTest,cTest$c.index)
  
}
######################################################################################################################################################


######################################################################################################################################################
####### 1.  penalized CC only Model ####################################################################################################################################
######################################################################################################################################################
predictiveModel_cox1 <- CoxModel$new()
predictiveModel_cox1$customTrain(RES,survObj)
fit<-predictiveModel_cox1$rawCaretModel()

cvResults1 <- crossValidatePredictiveCoxModel(RES, survObj, model=predictiveModel_cox1)

# by default, 5 fold cross validation
cindexTrain1<-c()
cindexTest1<-c()
for (i in 1:numFolds){
  cTrain<-concordance.index(x=cvResults1$trainPredictions[[i]], surv.time=cvResults1$trainObservations[[i]][,"time"], surv.event=cvResults1$trainObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTrain1 <- c(cindexTrain1,cTrain$c.index)
  cTest<-concordance.index(x=cvResults1$testPredictions[[i]], surv.time=cvResults1$testObservations[[i]][,"time"], surv.event=cvResults1$testObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTest1 <- c(cindexTest1,cTest$c.index)
}
######################################################################################################################################################


######################################################################################################################################################
####### 1-1.  which metric should be used in the penalized CC only Model #############################################################################
####### If any metric is proposed or given, then I can implement it with following unit test   #######################################################
######################################################################################################################################################
coeff1<-list()
coeff2<-list()
for(k in 1:numFolds){
  coeff1[[k]]<-findOptimalCoxModel(RES,survObj,model=myCoxModel$new(), Grid = createENetTuneGrid(alphas=1),numFolds = 10)
  coeff2[[k]]<-findOptimalCoxParam_DIndex(RES,survObj,model=myCoxModel$new(), Grid = createENetTuneGrid(alphas=1),numFolds = 10)
}


# quantile mean computation
IQM <- function(vec,na.rm = TRUE){
  a<-quantile(vec,na.rm = TRUE)
  b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
  return(mean(vec[b],na.rm = TRUE))
}

# To make a metric vector with all lambda 
INDEX <- function(coeff){
  c1<-list()
  c2<-list()
  # each lambda
  for (k in 1:length(coeff[[1]])){
    C1<-c()
    C2<-c()
  
    # cross validate
    for (kk in 1:length(coeff)){
      C1<-c(C1,coeff[[kk]][[k]]$CITrain)
      C2<-c(C2,coeff[[kk]][[k]]$CITest)
    }
    c1[[k]]<-C1
    c2[[k]]<-C2
  }
  
  # combine 
  CC1<-c()
  CC2<-c()
  for (k in 1:length(c1)){
    CC1<-c(CC1,IQM(c1[[k]],na.rm = TRUE))
    CC2<-c(CC2,IQM(c2[[k]],na.rm = TRUE))
  }
  return(list(trainIndex = CC1,testIndex = CC2))
}


CINDEX <- INDEX(coeff1)
DINDEX <- INDEX(coeff2)

par(mfrow = c(2,2))
plot(CINDEX$trainIndex,col = "blue",main = "C.index with Train")
plot(CINDEX$testIndex,col = "blue",main = "C.index with Test")

plot(DINDEX$trainIndex,col = "red",main = "D.index with Train")
plot(DINDEX$testIndex,col = "red",main = "D.index with Test")

