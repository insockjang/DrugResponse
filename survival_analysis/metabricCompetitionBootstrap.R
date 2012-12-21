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


nboots <- 100
x <- 1:(round(2*dim(FEA)[2] /3))    
theta <- function(x){sample(x,length(x),replace=TRUE)} 
results <- bootstrap(x,nboots,theta)

myCoxFunction <-function(trainInput,trainOutput,testInput,testOutput){
  
  # outputTrain/outputTest should be survObj  
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

Result <-list()
for(k in 1:10){

  training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))
  subResult <-list()
  for(kk in 1:nboots){
    Training <- training[results$thetastar[,kk]]
    survObjTrain <- Surv(as.numeric(dataSets_expr_clinical$responseData["survYears",Training]), as.numeric(dataSets_expr_clinical$responseData["survDeath",Training]))
    survObjTest <- Surv(as.numeric(dataSets_expr_clinical$responseData["survYears",-training]), as.numeric(dataSets_expr_clinical$responseData["survDeath",-training]))

    #FEA_training <- FEA[,Training]
    #FEA_testing <- FEA[,-training]
    trainInput1 <- t(RES[,Training])
    testInput1 <- t(RES[,-training])

  
    subResult[[kk]] <-myCoxFunction(trainInput1,survObjTrain,testInput1,survObjTest)
  }
  Result[[k]]<-subResult
}

############## you can select following training sets
# Model 2 : MF only
trainInput2 <- t(FEA_training) 
testInput2 <- t(FEA_testing) 

# Model 3 : MF+CC
trainInput3 <- t(rbind(FEA_training,RES_training)) 
testInput3 <- t(rbind(FEA_testing,RES_testing)) 
############################################################

################ simple analysis

COEF <-c()
for (k in 1:10){
  for (kk in 1:nboots){
    COEF <- cbind(COEF,as.numeric(Result[[k]][[kk]]$coefficients))
  }
}
apply(COEF !=0,1,sum)

cindexTrain <-c()
cindexTest <-c()
for (k in 1:10){
  for (kk in 1:nboots){
    cindexTrain <- c(cindexTrain,Result[[k]][[kk]]$cindexTrain$c.index)
    cindexTest <- c(cindexTest,Result[[k]][[kk]]$cindexTest$c.index)
  }
}

boxplot(cbind(cindexTrain,cindexTest))

