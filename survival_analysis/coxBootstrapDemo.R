library(predictiveModeling)
library(glmnet)
library(survival)
library(survcomp)
library(randomSurvivalForest)
library(synapseClient)

synapseLogin()

## competition : expression data
exprSet <- loadEntity(138993)

## competition : clinical covariate data
clinSet <- loadEntity(138996)

FEA <- exprs(exprSet$objects$eSet1)
RES <- exprs(clinSet$objects$eSet2)

## filter NA out from clinical datasets
a<-1:dim(RES)[2]
for(i in 1:dim(RES)[1]){
  a<-intersect(a,which(is.na(RES[i,])==0))
}
RES <- RES[,a]
FEA <- FEA[,a]


FEA<-t(scale(t(FEA)))

## how to compute
CONCORD <- function(X, Y, Z) { 
  tt <- concordance.index(x=X, surv.time=Y, surv.event=Z, method="noether", na.rm=TRUE); 
  return(   c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower,"upper"=tt$upper)) 
}


C_train1<-list()
C_test1<-list()
C_train2<-list()
C_test2<-list()
C_train3<-list()
C_test3<-list()

fit_cox1<-list()
fit_cox2<-list()
fit_cox3<-list()

predTest_cox1<-list()
predTest_cox2<-list()
predTest_cox3<-list()

predTrain_cox1<-list()
predTrain_cox2<-list()
predTrain_cox3<-list()

# Split data for training and testing
# For demo, I split the data 2/3 and 1/3, training and testing respectively.
for (kk in 1:10){
  training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))
  
  
  # after match sample annotation, we can split them into training and testing set: N fold cross validation
  nfolds <- 3
  nboots <- 100
  foldRatio <- (nfolds-1)/nfolds
  
  x <- 1:length(training)    
  
  theta <- function(x){sample(x,length(x),replace=TRUE)} 
  
  results <- bootstrap(x,nboots,theta)     
  
  bootSample<-results$thetastar
  
  fit_cox1[[kk]]<-list()
  fit_cox2[[kk]]<-list()
  fit_cox3[[kk]]<-list()
  
  predTest_cox1[[kk]]<-list()
  predTest_cox2[[kk]]<-list()
  predTest_cox3[[kk]]<-list()
  
  predTrain_cox1[[kk]]<-list()
  predTrain_cox2[[kk]]<-list()
  predTrain_cox3[[kk]]<-list()
  
  C_tr1 <- c()
  C_te1 <- c()
  C_tr2 <- c()
  C_te2 <- c()
  C_tr3 <- c()
  C_te3 <- c()
  
  for (k in 1:nboots){
    FEA_training <- FEA[,training[bootSample[,k]]]
    FEA_testing <- FEA[,-training]
    
    RES_training <- RES[,training[bootSample[,k]]]
    RES_testing <- RES[,-training]
    
    survObj<-Surv(time=RES_training[8,],event = RES_training[9,])
    survObj_test<-Surv(time=RES_testing[8,],event = RES_testing[9,])
    
    
    ###### lasso set alpha =1
    #### modified cox model $1
    #### CC + MF
    fit_cox1[[kk]][[k]] <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:7),])),survObj,family = "cox",alpha =1,nfolds=3)
    predTest_cox1[[kk]][[k]] <- predict(fit_cox1[[kk]][[k]],t(rbind(FEA_testing,RES_testing[c(1:7),])),s="lambda.min")
    predTrain_cox1[[kk]][[k]] <- predict(fit_cox1[[kk]][[k]],t(rbind(FEA_training,RES_training[c(1:7),])),s="lambda.min")
    C_tr1<-rbind(C_tr1,CONCORD(predTrain_cox1[[kk]][[k]],RES_training[8,],RES_training[9,]))
    C_te1<-rbind(C_te1,CONCORD(predTest_cox1[[kk]][[k]],RES_testing[8,],RES_testing[9,]))
    
    #### MF only : modified cox model $3
    fit_cox3[[kk]][[k]] <-cv.glmnet(t(FEA_training),survObj,family = "cox",alpha =1,nfolds=3)
    predTest_cox3[[kk]][[k]] <- predict(fit_cox3[[kk]][[k]],t(FEA_testing),s="lambda.min")
    predTrain_cox3[[kk]][[k]] <- predict(fit_cox3[[kk]][[k]],t(FEA_training),s="lambda.min")
    C_tr3<-rbind(C_tr3,CONCORD(predTrain_cox3[[kk]][[k]],RES_training[8,],RES_training[9,]))
    C_te3<-rbind(C_te3,CONCORD(predTest_cox3[[kk]][[k]],RES_testing[8,],RES_testing[9,]))
    
    #### CC only : modified cox model $2 
    # This approach means that each clinical covariate has its supplementary molecular features (selected by "feature selection by glmnet")
    # Then I combine all features and run cox model
    fit_cox2[[kk]][[k]] <-cv.glmnet(t(RES_training[c(1:7),]),survObj,family = "cox",alpha =1,nfolds=3)
    predTest_cox2[[kk]][[k]] <- predict(fit_cox2[[kk]][[k]],t(RES_testing[c(1:7),]),s="lambda.min")
    predTrain_cox2[[kk]][[k]] <- predict(fit_cox2[[kk]][[k]],t(RES_training[c(1:7),]),s="lambda.min")
    C_tr2<-rbind(C_tr2,CONCORD(predTrain_cox2[[kk]][[k]],RES_training[8,],RES_training[9,]))
    C_te2<-rbind(C_te2,CONCORD(predTest_cox2[[kk]][[k]],RES_testing[8,],RES_testing[9,]))
    
    print(k)
  }
  
  C_train1[[kk]]<-C_tr1
  C_test1[[kk]]<-C_te1
  
  C_train2[[kk]]<-C_tr2
  C_test2[[kk]]<-C_te2
  
  C_train3[[kk]]<-C_tr3
  C_test3[[kk]]<-C_te3
  
}


# check what kind of covariates and molecular features are significantly selected
coefMatrix1 <-c()

for (k in 1:10){
  coefMatrix <- matrix(0,nrow=(dim(FEA_training)[1] + dim(RES_training[c(1:7),])[1]),ncol=nboots)
  
  for (i in 1:nboots){
    coefMatrix[which(coef(fit_cox1[[k]][[i]])!=0),i]<-1
  }
  coefMatrix1 <- cbind(coefMatrix1,coefMatrix)
}
apply(coefMatrix1,1,sum)

# check what kind of covariates are significantly selected
coefMatrix2 <-c()

for (k in 1:10){
  coefMatrix <- matrix(0,nrow=dim(RES_training[c(1:7),])[1],ncol=nboots)
  
  for (i in 1:nboots){
    coefMatrix[which(coef(fit_cox2[[k]][[i]])!=0),i]<-1
  }
  coefMatrix2 <- cbind(coefMatrix2,coefMatrix)
}
apply(coefMatrix2,1,sum)



# check what kind of molecular features are significantly selected
coefMatrix3 <-c()

for (k in 1:10){
  coefMatrix <- matrix(0,nrow=dim(FEA_training)[1],ncol=nboots)
  
  for (i in 1:nboots){
    coefMatrix[which(coef(fit_cox3[[k]][[i]])!=0),i]<-1
  }
  coefMatrix3 <- cbind(coefMatrix3,coefMatrix)
}
apply(coefMatrix3,1,sum)



## plot cindex
cindex1_tr <-c()
cindex1_te <-c()

for (k in 1:10){  
  cindex1_tr<-c(cindex1_tr,C_train[[k]][,1])
  cindex1_te<-c(cindex1_te,C_test[[k]][,1])
}

boxplot(cbind(cindex1_tr,cindex1_te))

## plot cindex
cindex2_tr <-c()
cindex2_te <-c()

for (k in 1:10){  
  cindex2_tr<-c(cindex2_tr,C_train2[[k]][,1])
  cindex2_te<-c(cindex2_te,C_test2[[k]][,1])
}

boxplot(cbind(cindex2_tr,cindex2_te))

## plot cindex
cindex3_tr <-c()
cindex3_te <-c()

for (k in 1:10){  
  cindex3_tr<-c(cindex3_tr,C_train[[k]][,1])
  cindex3_te<-c(cindex3_te,C_test[[k]][,1])
}

boxplot(cbind(cindex3_tr,cindex3_te))