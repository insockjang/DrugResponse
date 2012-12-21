library(synapseClient)
synapseLogin()

# This code is implemented for curating data which has a features with gene symbol and non replicated samples with averaged replicates
# uploaded in Synapse with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
# specify your working directory where raw data exist
# Here I selected the Sanger dataset for curating
workDirectory <-  "/gluster/external-data/DAT_023__GSKCellLines/2010-12/Expression/Raw"
setwd(workDirectory)

library(affy)
library(preprocessCore)
# depending on CDF of given microarray, we have to use corresponding CDF info.
# Raw data importing procedure
Data <- ReadAffy()
# you can easily check the CDF by just typing
Data
# okay, GSK Expr microarray is hgu133plus2
library(hgu133plus2.db)
library(gdata)

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
x <- hgu133plus2SYMBOL
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
sumExp <-c()
for (i in 1:length(geneSymbolName)){
  probesPerGene<-probeGeneMap[which(!is.na(match(probeGeneMap[,2],geneSymbolName[i]))),1]
  set<-c()
  for(j in 1:length(probesPerGene)){
    set<-union(set,which(!is.na(match(probeName,probesPerGene[j]))))    
  }
  exp<-normPmDataBgRma[set,]
  sumExp <- rbind(sumExp,apply(t(exp)-medpolish(t(exp),trace.iter = FALSE)$residuals,1,mean))
}
rownames(sumExp) <- geneSymbolName

#########################################
## Read Sample Annotation for processing replicates
## Replace CEL file name to sampleName given annotation file
#########################################
AnnotData <- read.delim2("GSK_RNA.sdrf",header = TRUE, skip =8)

# Replace special characters "-" to "_"
annotSampleName <-toupper(AnnotData$Characteristics.Cell.Line.Name.)
while(length(grep(" ",annotSampleName)) !=0){
  annotSampleName<-sub(" ","",annotSampleName)
}
while(length(grep("-",annotSampleName)) !=0){
  annotSampleName<-sub("-","",annotSampleName)
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
    finalData<-cbind(finalData,sumExp[,a])  
  }
  else{
    finalData<-cbind(finalData,apply(sumExp[,a],1,gm.mean))
  }
}

colnames(finalData)<-uniqueSampleName
rownames(finalData)<-rownames(sumExp)


eSet <- new("ExpressionSet",exprs = finalData, annotation = annotation(Data))


probeIDLayer <- Layer(list(name = "R_Expression_GSK_RMA_GeneSymbol", type = "E", parentId = "114505", status="db"))
testLayer <- addObject(probeIDLayer, eSet)
testLayer <- storeEntity(testLayer)



