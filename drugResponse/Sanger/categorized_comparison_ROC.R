# DEMO naive Bayes Classification
library(predictiveModeling)
library(class)
library(e1071)
library(biclust)
library(clinfun)
library(randomForest)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")







source("/home/ijang/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")


###################################################
#### Load Sanger Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "210937"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "266141"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_exprLayer <- "210931" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "220680" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug


roc<-function(A){
  a<-matrix(0,nrow=2,ncol=2)
  for(kkk in 1:nrow(A)){
    if(A[kkk,1]==1 && A[kkk,2]==1) a[1,1]=a[1,1]+1
    if(A[kkk,1]==1 && A[kkk,2]==0) a[1,2]=a[1,2]+1
    if(A[kkk,1]==0 && A[kkk,2]==1) a[2,1]=a[2,1]+1
    if(A[kkk,1]==0 && A[kkk,2]==0) a[2,2]=a[2,2]+1    
  }
  
  TP<-a[1,1]
  FP<-a[1,2]
  FN<-a[2,1]
  TN<-a[2,2]
  TPR <-TP/(TP+FN)
  FPR <-FP/(FP+TN)  
  b<-c(TPR,FPR)
  names(b)<-c("TPR","FPR")
  return(b)
}

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

for(k1 in 1:nrow(adf_drug)){
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,1,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  threshold<-seq(min(filteredResponseDataScaled),max(filteredResponseDataScaled),len =100)
  
  threshold1<-seq(min(filteredResponseDataScaled),max(filteredResponseDataScaled),by =0.5)
  RF<-c()
  ##################################################################################### categorical results
      
  rfModel<-crossValidatePredictiveModelLogistic(filteredFeatureDataScaled,filteredResponseDataScaled,model = myRandomForestModel$new(),numFolds=5)    
  nbModel<-crossValidatePredictiveModelLogistic(feature,response,model = myNaiveBayesModel$new(),numFolds=5)  
    
  
  ######################################################################################### continuous 
  filename <-paste("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/ENet2/cvDrug_",k1,".Rdata",sep="")
  load(filename)
  
  pred<-c()
  obs <-c()
  for(kk in 1:length(resultsScale)){
    pred<-c(pred,resultsScale[[kk]]$testPredictions)
    obs<-c(obs,resultsScale[[kk]]$testObservations)
  }
  
  
  ##################################################################################### binarized results
  cont<-c()
  #threshold<-seq(min(min(obs),min(pred)),max(max(pred),max(obs)),len =100)  
  for(k in 21:99){
    OBS <- (obs<threshold[k])
    PRED <- (pred<threshold[k])    
    cont<-rbind(cont,roc(cbind(PRED,OBS)))
  }

  cont<-list()
  #  threshold<-seq(min(obs),max(obs),len =100)  
  for(k in 1:length(threshold)-2){
    OBS <- binarize(obs>=threshold[70])
    plot(roc.curve(pred,OBS))    
  }
  

   
  
  plot(RF[,2],RF[,1],col="red")
  points(NB[,2],NB[,1],col="blue")  
  points(cont[,2],cont[,1],col="green")
  abline(0,1)
  
}

plot(cont[,2],cont[,1],col="blue")