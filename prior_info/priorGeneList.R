library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
# source for priorList 
# 1. transcription factor list
# 2. Sanger Cancer Gene Census list

setwd("/Volumes/ijang/PredictiveModel/Regulator/")
TF<-read.table("tf-homo-current.txt")
entrezTF<-as.numeric(as.character(TF[,1]))
geneSymbolTF<-as.character(TF[,2])

M<-read.delim("sigenzymolome_cellular_feb10.txt")
entrezM<-as.numeric(as.character(M[,1]))
geneSymbolM<-as.character(M[,2])



setwd("/Volumes/ijang/PredictiveModel/CancerGeneCensus/")
library(gdata)
CGC<-read.xls("cancer_gene_census.xls")
cgc<- as.character(CGC$Symbol)

entrezCgc<- as.numeric(as.character(CGC$GeneID))

priorGeneList<-list(cancerGeneCensus = cgc,modulator = geneSymbolM, transcriptionFactor = geneSymbolTF)


dataset <- Dataset(list(
  name="prior knowledge gene list",
  parentId=170143
  ))
dataset <- createEntity(dataset)

#onWeb(dataset)

probeIDLayer <- Layer(list(name = "GeneSymbolList", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, priorGeneList)
testLayer <- storeEntity(testLayer)
