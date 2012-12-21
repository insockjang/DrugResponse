### DEMO Lasso
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")

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


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)

rankTransform <- function(FEA){
  FEA <- (apply(FEA,2,rank)-1)/(nrow(FEA)-1) - 0.5
}

Results<-list()

for(kk in 1:ncol(dataSets_ccle$responseData)){  
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE])
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  # rank transform
  filteredFeatureDataRank <- rankTransform(filteredFeatureData)
  filteredResponseDataRank <- (rank(filteredResponseData)-1)/(length(filteredResponseData)-1) - 0.5
  
  # 5 fold cross validation 
  alphas <- 1
  lambdas <- exp(seq(-5, 2, length = 100))
  
  
  resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
  resultsRank<-crossValidatePredictiveModel1(filteredFeatureDataRank, filteredResponseDataRank, model = myEnetModel$new(), alpha=alphas, lambda = lambdas, nfolds = 5)
  
  filename <- paste("NPLasso/cvDrug_",kk,".Rdata",sep ="")  
  save(resultsScale,resultsRank,lambdas,file = filename)
  Results[[kk]]<-list(resultsScale = resultsScale,resultsRank=resultsRank)
}


corPearsonScale<-c()
corSpearmanScale<-c()
corPearsonRank<-c()
corSpearmanRank<-c()

for(kk in 1:ncol(dataSets_ccle$responseData)){
  trPred <- foreach(k = 1:5) %do%{Results[[kk]]$resultsScale[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{Results[[kk]]$resultsScale[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{Results[[kk]]$resultsScale[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{Results[[kk]]$resultsScale[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  corPearsonScale<-rbind(corPearsonScale,c(cor(allTrPred,allTrObsr),cor(allTePred,allTeObsr)))
  corSpearmanScale<-rbind(corSpearmanScale,c(cor(allTrPred,allTrObsr,method = "spearman"),cor(allTePred,allTeObsr,method = "spearman")))
  
  trPred <- foreach(k = 1:5) %do%{Results[[kk]]$resultsRank[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{Results[[kk]]$resultsRank[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{Results[[kk]]$resultsRank[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{Results[[kk]]$resultsRank[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  corPearsonRank<-rbind(corPearsonRank,c(cor(allTrPred,allTrObsr),cor(allTePred,allTeObsr)))
  corSpearmanRank<-rbind(corSpearmanRank,c(cor(allTrPred,allTrObsr,method = "spearman"),cor(allTePred,allTeObsr,method = "spearman")))
  
}

filename1 <- "NPLasso/cv5performance.Rdata"

save(corPearsonScale,corSpearmanScale,corPearsonRank,corSpearmanRank,file = filename1)
