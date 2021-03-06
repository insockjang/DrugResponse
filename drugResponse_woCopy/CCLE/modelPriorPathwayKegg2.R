# simple comparison with overlap group lasso vs. lasso
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
# library(standGL)
source("~/COMPBIO/trunk/users/jang/R5/myNewSGLModel.R")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")
library(foreach)
# library(multicore)
# library(doMC)
# registerDoMC()
# load MsigDB at the symbol level
entity<-loadEntity(105350)
mSigDB_annotations <- loadEntity(105363)

DB<-entity$objects$MsigDB_symbolID
allPathways <- mSigDB_annotations$objects$C2$KEGG

# CCLE
# expressionSet
id_copyLayer <- "269019"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "269021"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_exprLayer <- "269056" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

# Because cnv data struggle with singualr problem when applying overlapping group lasso, I only use expressionSet and oncomapSet as my feature datasets
# ds_features_cn_mut_ccle <- createAggregateFeatureDataSet(list(expr = exprSet, copy = copySet,mut = oncomapSet))

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


for(kk in 1:ncol(dataSets_ccle$responseData)){  
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredFeatureData  <- t(unique(t(filteredFeatureData)))
  
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  DATA<-filteredFeatureDataScaled
  DATA1<-c()
  
  alphas = seq(0,1,len = 20)
  # When applying group lasso in overlaps, each group has its own index.
  #  
  index <-c()
  j<-0
  for (i in 1:length(allPathways)){  
    mSigDB_index1 <- which(DB$geneset.names == allPathways[[i]])
    curPathwayGenes1 <- DB$genesets[[mSigDB_index1]]
    a1<-paste(curPathwayGenes1,"_expr",sep="")
    a2<-paste(curPathwayGenes1,"_copy",sep="")
    a3<-paste(curPathwayGenes1,"_mut",sep="")
    A<-union(a1,union(a2,a3))
    AA<-colnames(DATA) %in% A
    pathwayGene <- DATA[,which(AA==1)]
    if (length(dim(pathwayGene)[1]) !=0){
      j<-j+1
      a<-dim(pathwayGene)[2]  
      b<-rep(j,a)
      index <-c(index,b)
      DATA1 <-cbind(DATA1,pathwayGene)
    }
  }
  data<-DATA[,setdiff(colnames(DATA),colnames(DATA1))]
  
  DATA2<-cbind(DATA1,data)
  index2<-seq(1,ncol(data))+max(index)
  index<-c(index,index2)
  # 5 fold cross validation 
  
  resultsScale<-crossValidatePredictiveModel1(DATA2, filteredResponseDataScaled, model = myNewSGLModel$new(), alpha=alphas, numFolds = 5,index = index,nfolds = 5)
  filename <- paste("priorPathwayKegg2/cvDrug_",kk,".Rdata",sep ="")  
  save(resultsScale,alphas,file = filename)
  
}
