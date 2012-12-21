library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)

# synapseLogin() ### not required if configured for automatic login

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

######################################################################################################################################################
#######  unpenalized CC only Model : Cox Regression #################################################################################################
######################################################################################################################################################
foldIndices <- createFolds(RES[,1],k=5)

RRR<-data.frame(res)
coxTrain<-list()
coxTest<-list()
for (k in 1:5){
  resTrain<-RRR[-foldIndices[[k]],]
  resTest<-RRR[foldIndices[[k]],]
  
  survObjTrain <-Surv(dataSets_expr_clinical$responseData[-foldIndices[[k]],"survYears"], dataSets_expr_clinical$responseData[-foldIndices[[k]],"survDeath"])
  survObjTest <-Surv(dataSets_expr_clinical$responseData[foldIndices[[k]],"survYears"], dataSets_expr_clinical$responseData[foldIndices[[k]],"survDeath"])
  fit<-coxph(survObjTrain ~ Site + ageDiagnosis + lymphnodes + grade + histology + tumorSizeCM + chemo + hormone + radiation + HER2 +ER + PR + ERPR + tripleNegative, resTrain)
  coxTrain[[k]] <- predict(fit)
  coxTest[[k]] <- predict(fit, resTest)
}
# by default, 5 fold cross validation
cindexTrain<-c()
cindexTest<-c()
for (i in 1:5){
  cTrain<-concordance.index(x=coxTrain[[i]], surv.time=survObj[-foldIndices[[i]],1], surv.event=survObj[-foldIndices[[i]],2], na.rm=TRUE, alpha= .05)
  cindexTrain <- c(cindexTrain,cTrain$c.index)
  cTest<-concordance.index(x=coxTest[[i]], surv.time=survObj[foldIndices[[i]],1], surv.event=survObj[foldIndices[[i]],2], na.rm=TRUE, alpha= .05)
  cindexTest <- c(cindexTest,cTest$c.index)
}
######################################################################################################################################################


######################################################################################################################################################
####### 1.  penalized CC only Model ####################################################################################################################################
######################################################################################################################################################
predictiveModel_cox1 <- CoxModel$new()
cvResults1 <- crossValidatePredictiveCoxModel(RES, survObj, model=predictiveModel_cox1)

# by default, 5 fold cross validation
cindexTrain1<-c()
cindexTest1<-c()
for (i in 1:5){
  cTrain<-concordance.index(x=cvResults1$trainPredictions[[i]], surv.time=cvResults1$trainObservations[[i]][,"time"], surv.event=cvResults1$trainObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTrain1 <- c(cindexTrain1,cTrain$c.index)
  cTest<-concordance.index(x=cvResults1$testPredictions[[i]], surv.time=cvResults1$testObservations[[i]][,"time"], surv.event=cvResults1$testObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTest1 <- c(cindexTest1,cTest$c.index)
}
######################################################################################################################################################


######################################################################################################################################################
####### 2. penalized MF only Model ####################################################################################################################################
######################################################################################################################################################
predictiveModel_cox1 <- CoxModel$new()
cvResults2 <- crossValidatePredictiveCoxModel(FEA, survObj, model=predictiveModel_cox1)

# by default, 5 fold cross validation
cindexTrain2<-c()
cindexTest2<-c()
for (i in 1:5){
  cTrain<-concordance.index(x=cvResults2$trainPredictions[[i]], surv.time=cvResults2$trainObservations[[i]][,"time"], surv.event=cvResults2$trainObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTrain2 <- c(cindexTrain2,cTrain$c.index)
  cTest<-concordance.index(x=cvResults2$testPredictions[[i]], surv.time=cvResults2$testObservations[[i]][,"time"], surv.event=cvResults2$testObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTest2 <- c(cindexTest2,cTest$c.index)
}
######################################################################################################################################################


######################################################################################################################################################
####### 3. penalized CC+MF Model ####################################################################################################################################
######################################################################################################################################################
predictiveModel_cox1 <- CoxModel$new()
cvResults3 <- crossValidatePredictiveCoxModel(cbind(FEA,RES), survObj, model=predictiveModel_cox1)

# by default, 5 fold cross validation
cindexTrain3<-c()
cindexTest3<-c()
for (i in 1:5){
  cTrain<-concordance.index(x=cvResults3$trainPredictions[[i]], surv.time=cvResults3$trainObservations[[i]][,"time"], surv.event=cvResults3$trainObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTrain3 <- c(cindexTrain3,cTrain$c.index)
  cTest<-concordance.index(x=cvResults3$testPredictions[[i]], surv.time=cvResults3$testObservations[[i]][,"time"], surv.event=cvResults3$testObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTest3 <- c(cindexTest3,cTest$c.index)
}
######################################################################################################################################################



######################################################################################################################################################
####### 4. unpenalized CC + penalized MF Model ####################################################################################################################################
######################################################################################################################################################
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myPenaltyCoxModel.R")
predictiveModel_cox2 <- myPenaltyCoxModel$new()
source("/home/ijang/COMPBIO/trunk/users/jang/survival_analysis/crossValidatePredictiveCoxPenaltyModel.R")
cvResults4 <- crossValidatePredictiveCoxPenaltyModel(cbind(FEA,RES), survObj, model=predictiveModel_cox2, penaltyFactor = c(rep(1,ncol(FEA)),rep(0,ncol(RES))))

# by default, 5 fold cross validation
cindexTrain4<-c()
cindexTest4<-c()
for (i in 1:5){
  cTrain<-concordance.index(x=cvResults4$trainPredictions[[i]], surv.time=cvResults4$trainObservations[[i]][,"time"], surv.event=cvResults4$trainObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTrain4 <- c(cindexTrain4,cTrain$c.index)
  cTest<-concordance.index(x=cvResults4$testPredictions[[i]], surv.time=cvResults4$testObservations[[i]][,"time"], surv.event=cvResults4$testObservations[[i]][,"status"], na.rm=TRUE, alpha= .05)
  cindexTest4 <- c(cindexTest4,cTest$c.index)
}
######################################################################################################################################################
######################################################################################################################################################



par(mfrow = c(1,2))
boxplot(cbind(cindexTrain,cindexTrain1,cindexTrain2,cindexTrain3,cindexTrain4),ylim = c(0.5,0.75),main = "Train")
boxplot(cbind(cindexTest,cindexTest1,cindexTest2,cindexTest3,cindexTest4),ylim = c(0.5,0.75),main = "Test")
dev.off()


Tr <-cbind(cindexTrain,cindexTrain1,cindexTrain2,cindexTrain3,cindexTrain4)

Te <- cbind(cindexTest,cindexTest1,cindexTest2,cindexTest3,cindexTest4)