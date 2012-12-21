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

featureData <- createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))

dataSets_expr_clinical <- createFeatureAndResponseDataList(featureData, clinicalData, filterNas = TRUE)

FEA <- t(scale(t(dataSets_expr_clinical$featureData)))

RES <- dataSets_expr_clinical$responseData

res<-matrix(0,nrow = 14,ncol = dim(RES)[2])
rownames(res)<-rownames(RES)[1:14]
colnames(res)<-colnames(RES)

res["Site",] <- as.numeric(RES["Site",])
res["ageDiagnosis",]<- as.numeric(RES["ageDiagnosis",])
res["lymphnodes",] <- ifelse(RES["lymphnodes",] == "pos",1,0)
res["grade",] <- as.numeric(RES["grade",])
res["tumorSizeCM",] <- as.numeric(RES["tumorSizeCM",])
res["chemo",] <- as.numeric(RES["chemo",])
res["hormone",] <- as.numeric(RES["hormone",])
res["radiation",] <- as.numeric(RES["radiation",])
res["HER2",] <- ifelse(RES["HER2",]=="pos", 1,0)
res["ER",] <- ifelse(RES["ER",]=="pos", 1,0)
res["PR",] <- ifelse(RES["PR",]=="pos", 1,0)
res["ERPR",] <- ifelse(RES["ERPR",]=="pos", 1,0)
res["tripleNegative",] <- ifelse(RES["tripleNegative",]=="pos", 1,0)
res["histology",] <- ifelse(RES["histology",]== "InfilitratingLobular",1,
                            ifelse(RES["histology",]== "InfiltratingDuctal",2,
                                   ifelse(RES["histology",]== "Medullary",3,
                                          ifelse(RES["histology",]== "MixedHistology",4,5))))
                     
rm(RES)
RES<-res

## factor
#res <- list()
#res$Site <- factor(RES["Site",])
#res$ageDiagnosis<- factor(RES["ageDiagnosis",])
#res$lymphnodes <- factor(RES["lymphnodes",])
#res$grade <- factor(RES["grade",])
#res$tumorSizeCM <- factor(RES["tumorSizeCM",])
#res$chemo <- factor(RES["chemo",])
#res$hormone <- factor(RES["hormone",])
#res$radiation <- factor(RES["radiation",])
#res$HER2 <- factor(RES["HER2",])
#res$ER <- factor(RES["ER",])
#res$PR <- factor(RES["PR",])
#res$ERPR <- factor(RES["ERPR",])
#res$tripleNegative <- factor(RES["tripleNegative",])
#res$histology <- factor(RES["histology",])
#RES<-data.frame(res)


myCoxFunction <-function(trainInput,trainOutput,testInput,testOutput){
  
  ######### LASSO COX model  
  myCoxModel <- CoxModel$new()
  myCoxModel$train(trainInput1, trainOutput)
  coefficients <- coef(myCoxModel$rawCaretModel(),s="lambda.min")
  myCoxModel_predTrain <- myCoxModel$customPredict(trainInput1)
  myCoxModel_predTest <- myCoxModel$customPredict(testInput1)
  
  cindexTrain = concordance.index(myCoxModel_predTrain,trainOutput[,1],trainOutput[,2],method ="noether")
  cindexTest = concordance.index(myCoxModel_predTest,testOutput[,1],testOutput[,2],method ="noether")
  
  return(list(coefficients = coefficients, cindexTrain=cindexTrain,cindexTest=cindexTest))
}

# 100 times for building boxplot for concordance index of training and testing(validating)
# data are randomly split into 2/3 training and 1/3 testing
# the reason why I implement this way is to run BigR for google gcomputing
# for loop might be easily exchanged into "foreach" and "%dopar%"
Result <- list()
for (k in 1:100){

  training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))

  survObjTrain <- Surv(as.numeric(dataSets_expr_clinical$responseData["survYears",training]), as.numeric(dataSets_expr_clinical$responseData["survDeath",training]))
  survObjTest <- Surv(as.numeric(dataSets_expr_clinical$responseData["survYears",-training]), as.numeric(dataSets_expr_clinical$responseData["survDeath",-training]))

  #### here I run CC-only model 
  #FEA_training <- FEA[,training]
  #FEA_testing <- FEA[,-training]

  RES_training <- RES[,training]
  RES_testing <- RES[,-training]

  trainInput1 <- t(RES_training) 
  testInput1 <- t(RES_testing) 

  Result[[k]]<-myCoxFunction(trainInput1,survObjTrain,testInput1,survObjTest)
}

### capture concordance index
cindexTrain <-c()
cindexTest <-c()
for (k in 1:100){
   cindexTrain <- c(cindexTrain,Result[[k]]$cindexTrain$c.index)
   cindexTest <- c(cindexTest,Result[[k]]$cindexTest$c.index)
}

## show the concordance index
boxplot(cbind(cindexTrain,cindexTest))

###################################################################################
########## below part is for single run

######### LASSO COX model
# Model 1 : CC only
trainInput1 <- t(RES_training) 
testInput1 <- t(RES_testing) 

myCoxModel1 <- CoxModel$new()
myCoxModel1$train(trainInput1, survObjTrain)
myCoxModel1_predTrain <- myCoxModel1$customPredict(trainInput1)
myCoxModel1_predTest <- myCoxModel1$customPredict(testInput1)

# Model 2 : MF only
trainInput2 <- t(FEA_training) 
testInput2 <- t(FEA_testing) 

myCoxModel2 <- CoxModel$new()
myCoxModel2$train(trainInput2, survObjTrain)
myCoxModel2_predTrain <- myCoxModel2$customPredict(trainInput2)
myCoxModel2_predTest <- myCoxModel2$customPredict(testInput2)

# Model 3 : MF+CC
trainInput3 <- t(rbind(FEA_training,RES_training)) 
testInput3 <- t(rbind(FEA_testing,RES_testing)) 

myCoxModel3 <- CoxModel$new()
myCoxModel3$train(trainInput3, survObjTrain)
myCoxModel3_predTrain <- myCoxModel3$customPredict(trainInput3)
myCoxModel3_predTest <- myCoxModel3$customPredict(testInput3)


############### analysis step

c1<-concordance.index(myCoxModel1_predTrain,survObjTrain[,1],survObjTrain[,2],method ="noether")
c2<-concordance.index(myCoxModel2_predTrain,survObjTrain[,1],survObjTrain[,2],method ="noether")
c3<-concordance.index(myCoxModel3_predTrain,survObjTrain[,1],survObjTrain[,2],method ="noether")

cindex.comp(c1,c2)
cindex.comp(c1,c3)
cindex.comp(c2,c3)

tc1<-concordance.index(myCoxModel1_predTest,survObjTest[,1],survObjTest[,2],method ="noether")
tc2<-concordance.index(myCoxModel2_predTest,survObjTest[,1],survObjTest[,2],method ="noether")
tc3<-concordance.index(myCoxModel3_predTest,survObjTest[,1],survObjTest[,2],method ="noether")
cindex.comp(tc1,tc2)
cindex.comp(tc1,tc3)
cindex.comp(tc2,tc3)


# you need to use forestplot for comparing models

textTrain <- t(t(c("ClinicalCovariate Only","Molecular Feature Only","Molecular Feature and Clinical Covariate")))
m <-c(c1$c.index,c2$c.index,c3$c.index)
l <- c(c1$lower,c2$lower,c3$lower)
u <- c(c1$upper,c2$upper,c3$upper)
forestplot(textTrain,m,l,u,zero=0.5,clip=c(log(0.1),log(2.5)), boxsize=0.75,xlab ="Concordance Index",
           col=meta.colors(box="royalblue",line="darkblue", summary="royalblue", zero = "darkred"))

textTest <- t(t(c("ClinicalCovariate Only","Molecular Feature Only","Molecular Feature and Clinical Covariate")))
m <-c(tc1$c.index,tc2$c.index,tc3$c.index)
l <- c(tc1$lower,tc2$lower,tc3$lower)
u <- c(tc1$upper,tc2$upper,tc3$upper)
forestplot(textTest,m,l,u,zero=0.5,clip=c(log(0.1),log(2.5)), boxsize=0.75,xlab ="Concordance Index",
           col=meta.colors(box="royalblue",line="darkblue", summary="royalblue", zero = "darkred"))

