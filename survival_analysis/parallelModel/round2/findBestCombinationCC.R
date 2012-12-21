library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@") ### not required if configured for automatic login

idExpressionLayer <- "160776" ##"160644" 
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idCopyLayer <- "160778" ##"160646"
copyLayer <- loadEntity(idCopyLayer)
copyData <- copyLayer$objects[[1]]

idClinicalFeaturesLayer <- "160780" ##"139171"
clinicalFeaturesLayer <- loadEntity(idClinicalFeaturesLayer)
clinicalFeaturesData <- clinicalFeaturesLayer$objects[[1]]

idClinicalSurvLayer <- "160782"
clinicalSurvLayer <- loadEntity(idClinicalSurvLayer)
clinicalSurvData <- clinicalSurvLayer$objects[[1]]


#### prepare feature data for predictive modeling by transposing the matrix to have samples on the rows and features on the columns and scaling the columns

featureData <-t(createAggregateFeatureDataSet(list(expr = exprData, copy = copyData)))
featureData_filtered <- filterNasFromMatrix(dataMatrix=featureData, filterBy = "columns")
featureData <- scale(featureData_filtered)

dataSets_MF_clinical <- createFeatureAndResponseDataList(featureData, clinicalSurvData)
dataSets_CC_clinical <- createFeatureAndResponseDataList(pData(clinicalFeaturesData), clinicalSurvData)

RES<- dataSets_CC_clinical$featureData

# clinical covariate numeric mapping
res<-matrix(0,ncol = 14,nrow = dim(RES)[1])
rownames(res)<-rownames(RES)
colnames(res)<-colnames(RES)[1:14]

res[,"Site"] <- as.numeric(RES[,"Site"])
res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
res[,"grade"] <- as.numeric(RES[,"grade"])
res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                            ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                   ifelse(RES[,"histology"]== "Medullary",3,
                                          ifelse(RES[,"histology"]== "MixedHistology",4,5))))

res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
res[,"chemo"] <- as.numeric(RES[,"chemo"])
res[,"hormone"] <- as.numeric(RES[,"hormone"])
res[,"radiation"] <- as.numeric(RES[,"radiation"])
res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
res[,"tripleNegative"] <- as.numeric(RES[,"tripleNegative"])

RES<-res

survObj <- dataSets_CC_clinical$responseData


######################################################################################################################################################
#######  Find the best combinations of Clinical Covariates in Cox Regression #################################################################################################
######################################################################################################################################################
# The number of folds for cross validation

binaryTable <-function(data){
  featureNum <- ncol(data)
  a<-c(0,1)
  bin<-expand.grid(a,a,a,a,a,a,a,a,a,a,a,a,a,a) # This might be changed when you have different number of clinical covariates
  names(bin)<-colnames(data)
  return(bin)
}

bT<-binaryTable(RES)
dRES<-data.frame(res)

ConcordanceCombinations<-list()
for (j in 1:10){
  foldIndices <- createFolds(dRES[,1],k=5)
  
  ConcordanceCombination <-foreach(fold = foldIndices) %dopar%{
    CTrain_survival<-c()
    CTest_survival<-c()
    
    for(k in 1:nrow(bT)){    
      if(sum(bT[k,])<=1){
        CTrain_survival <-c(CTrain_survival,0.5)
        CTest_survival <-c(CTest_survival,0.5)
      }
      else{
        fit <- coxph(survObj[-fold]~.,dRES[-fold,which(bT[k,]!=0)])
        CTrain_survival<-c(CTrain_survival,survConcordance(survObj[-fold] ~ predict(fit))$concordance)
        CTest_survival<-c(CTest_survival,survConcordance(survObj[fold] ~ predict(fit,dRES[fold,which(bT[k,]!=0)]))$concordance)      
      }
    }
    return(list(CTrain_survival = CTrain_survival,
                CTest_survival = CTest_survival))
  }
  ConcordanceCombinations[[j]]<-ConcordanceCombination
}

ctrain_survival <-c()
ctest_survival <-c()
for(k1 in 1:10){
  for(k in 1:5){
    ctrain_survival<-cbind(ctrain_survival,ConcordanceCombinations[[k1]][[k]]$CTrain_survival)
    ctest_survival<-cbind(ctest_survival,ConcordanceCombinations[[k1]][[k]]$CTest_survival)
  }
}



IQM <- function(vec,na.rm = TRUE){
  a<-quantile(vec,na.rm = TRUE)
  b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
  return(mean(vec[b],na.rm = TRUE))
}

CTestIQM<-apply(ctest_survival,1,IQM)
CTestMean<-apply(ctest_survival,1,mean)

CTrainIQM<-apply(ctrain_survival,1,IQM)
CTrainMean<-apply(ctrain_survival,1,mean)


which.max(CTestIQM)
which.max(CTestMean)

s<- sort(CTestMean,decreasing = TRUE,index.return = TRUE)
s$x[1:10]
bT[s$ix[1:100],]

apply(bT[s$ix[1:10],],2,sum)
range(s$x[1:10])

bT[which.max(CTestMean),]

# find the best combination in binary table which show me all of possible combination
best_IQM_combination<-bT[which(CTestIQM==max(CTestIQM)),]
best_mean_combination<-bT[which(CTestMean==max(CTestMean)),]


plot_colors <- c("blue","red","forestgreen","cyan","magenta")
par(mfrow = c(1,2))
plot(CTestMean,col = plot_colors[1],ylim = c(0.45,0.7),pch = 1,lty =1, axe = T,cex = .6,ylab = "Concordance Index",main = "Testing with Mean",xlab = "Possible Combinations")
ii<-which(CTestMean == max(CTestMean))
text(ii,CTestMean[ii],as.character(ii),offset = -0.1,cex = 1,pos =1)
points(ii,CTestMean[ii],pch = 18,cex = 2)
bT[which(CTestMean == max(CTestMean)),]

plot(CTestIQM,col = plot_colors[2],ylim = c(0.45,0.7),pch = 1,lty =1, axe = T,cex = .6,ylab = "Concordance Index",main = "Testing with IQM",xlab = "Possible Combinations")
ii<-which(CTestIQM == max(CTestIQM))
points(ii,CTestIQM[ii],pch = 18,cex = 2)
text(ii,CTestIQM[ii],as.character(ii),offset = -0.1)
bT[which(CTestIQM == max(CTestIQM)),]

bestCC <-bT[which(CTestIQM == max(CTestIQM)),]