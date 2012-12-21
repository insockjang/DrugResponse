# Sample code MsigDB Curation with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
## load the synapse client and login
library(synapseClient)
library(GSA)
library(affy)
synapseLogin()

## set up a project
myName <- "MsigDB all inclusive"
MsigDB_symbolID <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/msigdb.v3.0.symbols.gmt")
for(k in 1:length(MsigDB_symbolID$genesets)){
  a<-MsigDB_symbolID$genesets[[k]]
  b<-strsplit(gsub(" ","",a),"///")
  cc<-c()
  for(k1 in 1:length(b)){
    cc<-union(cc,b[[k1]])
  }
  MsigDB_symbolID$genesets[[k]]<-cc
}

MsigDB_entrezID <- GSA.read.gmt("msigdb.v3.0.entrez.gmt")
MsigDB_probeID <- GSA.read.gmt("msigdb.v3.0.orig.gmt")

projName <- sprintf("%ss MsigDB Project %s", myName, as.character(gsub("-",".",Sys.Date())))

myProj <- createEntity(Project(list(name=projName)))

dsProp <- list(name = "MsigDB",
               
               description = "BroadInstitute MsigDB is used to test predictive features. Fisher Exact Test will be conducted.",
               
               parentId = propertyValue(myProj, "id"))

dsAnnot <- list()

dataset <- Dataset(dsProp)

dataset <- createEntity(dataset)


onWeb(dataset)




probeIDLayer <- Layer(list(name = "MsigDB_probeID", type = "E", parentId = "105191", status="db"))
testFiles <- "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_probeID.ROBJECT"
testLayer <- addObject(probeIDLayer, MsigDB_probeID)
testLayer <- storeEntity(testLayer)


symbolIDLayer <- Layer(list(name = "MsigDB_symbolID_new", type = "E", parentId = "105191", status="db"))
testFiles <- "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_symbolID.ROBJECT"
load(testFiles)
testLayer2 <- addObject(symbolIDLayer, MsigDB_symbolID)
testLayer2 <- storeEntity(testLayer2)


entrezIDLayer <- Layer(list(name = "MsigDB_entrezID", type = "E", parentId = "105191", status="db"))
testFiles <- "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_entrezID.ROBJECT"
load(testFiles)
testLayer3 <- addObject(entrezIDLayer, MsigDB_entrezID)
testLayer3 <- storeEntity(testLayer3)


annotLayer <- Layer(list(name = "MsigDB_annotation", type = "E", parentId = "105191", status="db"))
testFiles <- "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_annotation.ROBJECT"
Annotation <-c()
Annotation$C1<-C1
Annotation$C2<-C2
Annotation$C3<-C3
Annotation$C4<-C4
Annotation$C5<-C5
load(testFiles)
testLayer4 <- addObject(annotLayer, Annotation)
testLayer4 <- storeEntity(testLayer4)


# original GMT file formats
probeIDLayer <- Layer(list(name = "GMT_MsigDB probeID", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testFiles <- "/Users/insockjang/Downloads/msigdb.v3.0.orig.gmt"
testLayer <- addFile(probeIDLayer, testFiles)
testLayer <- storeEntity(testLayer)


symbolIDLayer <- Layer(list(name = "GMT_MsigDB symbolID", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testFiles <- "/Users/insockjang/Downloads/msigdb.v3.0.symbols.gmt"
testLayer2 <- addFile(symbolIDLayer, testFiles)
testLayer2 <- storeEntity(testLayer2)


entrezIDLayer <- Layer(list(name = "GMT_MsigDB entrezID", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testFiles <- "/Users/insockjang/Downloads/msigdb.v3.0.entrez.gmt"
testLayer3 <- addFile(entrezIDLayer, testFiles)
testLayer3 <- storeEntity(testLayer3)

