# This code is implemented for curating data which has a features with gene symbol and non replicated samples with averaged replicates
# uploaded in Synapse with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
# specify your working directory where raw data exist
# Here I selected the Sanger dataset for curating
workDirectory <- "/Volumes/ijang/TEMP/test_folder/"
setwd(workDirectory)

library(affy)
library(preprocessCore)
# depending on CDF of given microarray, we have to use corresponding CDF info.
# Raw data importing procedure
Data <- ReadAffy()
# you can easily check the CDF by just typing
Data
# okay, Sanger Expr microarray is ht-hgu133a
library(hgu133a.db)
library(gdata)
#library(hgu133plus2.db)


probeName<-probeNames(Data)
probeSetName <- geneNames(Data)
## backgound adjustment : RMA convolution
# you might want to MAS5.0 background.Then replace method parameter to : method = "mas"
DataBgRma <- bg.correct(Data, method = 'rma')


## Warning !!! hereafter, data should be log2 transformed
# log2 transformation
pmDataBgRma<-log2(probes(DataBgRma,"pm"))

## normalization : Quantile Normalization
normPmDataBgRma<-normalize.quantiles(pmDataBgRma)
rownames(normPmDataBgRma) <-rownames(pmDataBgRma)
colnames(normPmDataBgRma) <-colnames(pmDataBgRma)

# Symbol and probeset mapping ###########################
x <- hgu133aSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
probeGeneMap <-c()
for(i in 1:length(xx)){
  probeGeneMap <- rbind(probeGeneMap,cbind(names(xx[i]),xx[[i]]))
}

#####################################
geneSymbolName<-unique(probeGeneMap[,2])

## Summarization : median polish by Tukey
# map probes to probeSetID
sumExp <-c()
for (i in 1:length(probeSetName)){
  set<-which(!is.na(match(probeName,probeSetName[i])))
  exp<-normPmDataBgRma[set,]
  sumExp <- rbind(sumExp,apply(t(exp)-medpolish(t(exp),trace.iter = FALSE)$residuals,1,mean))
}
rownames(sumExp) <- probeSetName

cv <- function(data){
  sqrt(var(data))/mean(data)  
}


# map the probeSetID having the highest coefficient of variation to geneSymbolID
reExp <-c()
probeMap <-c()
for (i in 1:length(geneSymbolName)){
  set<-probeGeneMap[which(!is.na(match(probeGeneMap[,2],geneSymbolName[i]))),1]
  partData<-sumExp[which(!is.na(match(probeSetName,set))),]
  if(length(set)==1){
    reExp<-rbind(reExp,sumExp[which(!is.na(match(probeSetName,set))),])
    probeMap<-cbind(probeMap,set)        
  }
  else{    
    reExp<-rbind(reExp,sumExp[which(!is.na(match(probeSetName,names(which.max(apply(partData,1,cv)))))),])
    probeMap<-cbind(probeMap,names(which.max(apply(partData,1,cv))))    
  }
}
rownames(reExp) <- geneSymbolName

#########################################
## Read Sample Annotation for processing replicates
## Replace CEL file name to sampleName given annotation file
#########################################
AnnotData <- read.xls("Sanger_affy_n798_sample_info_published.xls")

# Replace special characters "-" to "_"
annotSampleName <-toupper(AnnotData$SampleName)
while(length(grep("-",annotSampleName)) !=0){
  annotSampleName<-sub("-","_",annotSampleName)
}
while(length(grep("_",annotSampleName)) !=0){
  annotSampleName<-sub("_","",annotSampleName)
}

# define function : geometric mean, not arithmetic mean : because of log2 transformed data 
gm.mean <-function(dataVec){
  n <- length(dataVec)
  prod(dataVec) ^(1/n)  
}

# to preserve non replicated samples 
uniqueSampleName <-unique(annotSampleName)
finalData<-c()
for (i in 1:length(uniqueSampleName)){
  a<-which(uniqueSampleName[i]==annotSampleName)
  if (length(a)==1){
    finalData<-cbind(finalData,reExp[,a])  
  }
  else{
    finalData<-cbind(finalData,apply(reExp[,a],1,gm.mean))
  }
}

colnames(finalData)<-uniqueSampleName
rownames(finalData)<-rownames(reExp)


eSet <- new("ExpressionSet",exprs = finalData, annotation = annotation(Data))


probeIDLayer <- Layer(list(name = "R_Expression_Sanger_RMA_HCoV_GeneSymbol", type = "E", parentId = "114509", status="db"))
testLayer <- addObject(probeIDLayer, eSet)
testLayer <- storeEntity(testLayer)


