library(predictiveModeling)
library(glmnet)
library(survival)
library(survcomp)
library(randomSurvivalForest)
library(synapseClient)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

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

# Split data for training and testing
# For demo, I split the data 2/3 and 1/3, training and testing respectively.
training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))


# after match sample annotation, we can split them into training and testing set: N fold cross validation
nfolds <- 3
nboots <- 100
foldRatio <- (nfolds-1)/nfolds

x <- 1:length(training)    

theta <- function(x){sample(x,length(x),replace=TRUE)} 

results <- bootstrap(x,nboots,theta)     

bootSample<-results$thetastar

res <- list(Age35 = RES[1,],T = RES[2,],N=RES[3,],ERPR = RES[4,],HER2gain = RES[5,],Histology = RES[6,],grade = RES[7,],Time = RES[8,],Death = RES[9,])

## how to compute
CONCORD <- function(X, Y, Z) { 
  tt <- concordance.index(x=X, surv.time=Y, surv.event=Z, method="noether", na.rm=TRUE); 
  return(   c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower,"upper"=tt$upper)) 
}

fit_cox2<-list()
predTest_cox2<-list()
predTrain_cox2<-list()
C_tr2<-c()
C_te2<-c()

for (k in 1:nboots){
  res_training <- list(Age35 = RES[1,training[bootSample[,k]]],
                       T = RES[2,training[bootSample[,k]]],
                       N=RES[3,training[bootSample[,k]]],
                       ERPR = RES[4,training[bootSample[,k]]],
                       HER2gain = RES[5,training[bootSample[,k]]],
                       Histology = RES[6,training[bootSample[,k]]],
                       grade = RES[7,training[bootSample[,k]]],
                       Time = RES[8,training[bootSample[,k]]],
                       Death = RES[9,training[bootSample[,k]]]
                       )
  
  res_testing <- list(Age35 = RES[1,-training],
              T = RES[2,-training],
              N = RES[3,-training],
              ERPR = RES[4,-training],
              HER2gain = RES[5,-training],
              Histology = RES[6,-training],
              grade = RES[7,-training],
              Time = RES[8,-training],
              Death = RES[9,-training]
              )
  
  survObj<-Surv(time=res_training$Time,event = res_training$Death)
  survObj_test<-Surv(time=res_testing$Time,event = res_testing$Death)
  
  #### CC only : cox model $2 
  # Then I combine all clinical covariates(risk factors) and run cox model
  #fit_cox2[[k]] <-coxph(survObj ~ Age35+T+N+ERPR+HER2gain+Histology + grade, res_training)
  #fit_cox2[[k]] <-coxph(survObj ~ ERPR+HER2gain, res_training)
  fit_cox2[[k]] <-coxph(survObj ~ N, res_training)
  
  predTrain_cox2[[k]]<-predict(fit_cox2[[k]])
  predTest_cox2[[k]]<-predict(fit_cox2[[k]],res_testing)
  
  C_tr2<-rbind(C_tr2,CONCORD(predTrain_cox2[[k]],res_training$Time,res_training$Death))
  C_te2<-rbind(C_te2,CONCORD(predTest_cox2[[k]],res_testing$Time,res_testing$Death))
  
  
}

boxplot(cbind(C_tr2[,1],C_te2[,1]))

