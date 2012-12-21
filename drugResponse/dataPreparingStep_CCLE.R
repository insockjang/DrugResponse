library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
library(affy)

id_copyLayer <- "48339"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$copySet

id_oncomapLayer <- "48341"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$oncomapSet

id_exprLayer <- "48344" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$exprSet


cellLineName<-function(eSet){
  A<-colnames(exprs(eSet))
  B<-strsplit(as.character(A),"_")
  
  BB<-c()
  for(k in 1:length(B)){
    BB<-c(BB,B[[k]][1])
  }
  return(BB)
}

exprNames<-cellLineName(eSet_expr)
copyNames<-cellLineName(eSet_copy)
oncoNames<-cellLineName(eSet_oncomap)

exprNames[which(duplicated(exprNames)==1)] #778  #TT_OESOPHAGUS
copyNames[which(duplicated(copyNames)==1)] #849  #TT_OESOPHAGUS
oncoNames[which(duplicated(oncoNames)==1)] #756  #TT_OESOPHAGUS

colnames(exprs(eSet_expr))[grep("TT_",colnames(exprs(eSet_expr)))]
colnames(exprs(eSet_copy))[grep("TT_",colnames(exprs(eSet_copy)))]
colnames(exprs(eSet_oncomap))[grep("TT_",colnames(exprs(eSet_oncomap)))]

eSet_expr<-eSet_expr[,-778]
eSet_copy<-eSet_copy[,-849]
eSet_oncomap<-eSet_oncomap[,-756]

sampleNames(eSet_expr)<-toupper(exprNames[-778])
sampleNames(eSet_copy)<-toupper(copyNames[-849])
sampleNames(eSet_oncomap)<-toupper(oncoNames[-756])

##################################################################### Start of New table of drug
file = "/home/ijang/PredictiveModel/CCLE/CCLE_NP24.2009_Drug_data_2012.02.20.csv"

ccle<-read.csv(file, header = TRUE, sep = ",", quote="\"", dec=".",fill = TRUE)

A<-as.character(ccle$CCLE.Cell.Line.Name)
B<-strsplit(as.character(A),"_")

BB<-c()
for(k in 1:length(B)){
  BB<-c(BB,B[[k]][1])
}

BBB<-unique(BB)

CCC<-table(ccle$Compound)

CC<-as.character(ccle$Compound)
CCC<-unique(CC)

drugMatrix<-matrix(NA,nrow = length(BBB),ncol = length(CCC))
for(k1 in 1:length(BBB)){
  bb<-which(BB==BBB[k1])
  for(k2 in 1:length(CCC)){
    cc<-which(CC==CCC[k2])
    dd<-intersect(bb,cc)
    if(length(dd)==0){
      next
    }
    else{
      drugMatrix[k1,k2]<-as.numeric(ccle$ActArea[dd])
    }
  }
}

rownames(drugMatrix)<-make.names(toupper(BBB))
colnames(drugMatrix)<-(CCC)

adf_drug<-AnnotatedDataFrame(data=as.data.frame(drugMatrix))
##################################################################### End of New table of drug

myName <- "In Sock Jang"
project <- Project(list(
  name=paste("Cancer Cell Line Project - ", myName)
  ))
project <- createEntity(project)

onWeb(project)
# 
# analysis <- Analysis(list(
#   name="Elastic net versus custom methods",
#   description="Several Machine Learning methods run upon CCLE Data with Sanger Drug Response",
#   parentId=propertyValue(project, "id")
#   ))
# analysis <- createEntity(analysis)

dataset <- Dataset(list(
  name="CCLE",
  parentId=propertyValue(project, "id")
  ))
dataset <- createEntity(dataset)

#onWeb(dataset)

probeIDLayer <- Layer(list(name = "R_Expr", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_expr)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Cnv", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_copy)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Oncomap", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_oncomap)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Drug", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, adf_drug)
testLayer <- storeEntity(testLayer)
