library(predictiveModeling)
library(glmnet)
library(survival)
library(survcomp)
library(synapseClient)

synapseLogin()

## competition : expression data
exprSet <- loadEntity(138993)

## competition : clinical covariate dataq
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

FEA <- t(scale(t(FEA)))
# Split data for training and testing
# For demo, I split the data 2/3 and 1/3, training and testing respectively.
training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))

# 
FEA_training <- FEA[1:100,training]
FEA_testing <- FEA[1:100,-training]

RES_training <- RES[,training]
RES_testing <- RES[,-training]

survObj<-Surv(time=RES_training[8,],event = RES_training[9,])

######### LASSO COX model
# Model 1 : CC only
trainInput1 <- t(RES_training[1:7,]) 
testInput1 <- t(RES_testing[1:7,]) 

myCoxModel1 <- CoxModel$new()
myCoxModel1$train(trainInput1, survObj)
myCoxModel1_predTrain <- myCoxModel1$customPredict(trainInput1)
myCoxModel1_predTest <- myCoxModel1$customPredict(testInput1)

# Model 2 : MF only
trainInput2 <- t(FEA_training) 
testInput2 <- t(FEA_testing) 

myCoxModel2 <- CoxModel$new()
myCoxModel2$train(trainInput2, survObj)
myCoxModel2_predTrain <- myCoxModel2$customPredict(trainInput2)
myCoxModel2_predTest <- myCoxModel2$customPredict(testInput2)

# Model 3 : MF+CC
trainInput3 <- t(rbind(FEA_training,RES_training[1:7,])) 
testInput3 <- t(rbind(FEA_testing,RES_testing[1:7,])) 

myCoxModel3 <- CoxModel$new()
myCoxModel3$train(trainInput3, survObj)
myCoxModel3_predTrain <- myCoxModel3$customPredict(trainInput3)
myCoxModel3_predTest <- myCoxModel3$customPredict(testInput3)



##### PCR 
# pc_train <- svd(FEA_training - rowMeans(FEA_training))
# pc_test <- pc$U %*% pc$d * (FEA_testing - rowMeans(FEA_testing))

# cox.fit <- coxph(survObj ~ pc$v[,1:3])
# predict(cox.fit, pc_test[,1:3]])


### check/compare the performance among models
c1<-concordance.index(myCoxModel1_predTrain,RES_training[8,],RES_training[9,],method ="noether")
c2<-concordance.index(myCoxModel2_predTrain,RES_training[8,],RES_training[9,],method ="noether")
c3<-concordance.index(myCoxModel3_predTrain,RES_training[8,],RES_training[9,],method ="noether")

cindex.comp(c1,c2)
cindex.comp(c1,c3)
cindex.comp(c2,c3)

tc1<-concordance.index(myCoxModel1_predTest,RES_testing[8,],RES_testing[9,],method ="noether")
tc2<-concordance.index(myCoxModel2_predTest,RES_testing[8,],RES_testing[9,],method ="noether")
tc3<-concordance.index(myCoxModel3_predTest,RES_testing[8,],RES_testing[9,],method ="noether")
cindex.comp(tc1,tc2)
cindex.comp(tc1,tc3)
cindex.comp(tc2,tc3)


# you need to use forestplot for comparing models
#par(mfrow = c(2,1))
textTrain <- t(t(c("ClinicalCovariate Only","Molecular Feature Only","Molecular Feature and Clinical Covariate")))
m <-c(c1$c.index,c2$c.index,c3$c.index)
l <- c(c1$lower,c2$lower,c3$lower)
u <- c(c1$upper,c2$upper,c3$upper)
forestplot(textTrain,m,l,u,zero=0,clip=c(log(0.1),log(2.5)), boxsize=0.75,
           col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))

textTest <- t(t(c("ClinicalCovariate Only","Molecular Feature Only","Molecular Feature and Clinical Covariate")))
m <-c(tc1$c.index,tc2$c.index,tc3$c.index)
l <- c(tc1$lower,tc2$lower,tc3$lower)
u <- c(tc1$upper,tc2$upper,tc3$upper)
forestplot(textTest,m,l,u,zero=0,clip=c(log(0.1),log(2.5)), boxsize=0.75,
           col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))
#dev.off()

