library(foreach)
library(synapseClient)
library(predictiveModeling)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  

allPathways <- mSigDB_annotations$objects$C2$BIOCARTA

# preparing for FET : make total set
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
allGenes <-c()
for (i in 1:length(geneAllSetList)){
  allGenes<-union(allGenes,geneAllSetList[[i]])
}

Folder <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/NP/Ridge2/"

#AllResults<-list()

for(kk in 1:24){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(Folder,"fsDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  
  MAT<-as.matrix(resultsRidge)
  NAME <- rownames(MAT)
  MAT1<-MAT[-1]
  names(MAT1)<-NAME[-1]
  referenceSet<-sort(abs(MAT1), decreasing =T, index.return =T)
  
  
  testSet = names(referenceSet$x[1:500])
  testSet1 = strsplit(testSet,"_")
  testSet2<-c()
  for(k in 1:length(testSet1)){
    testSet2<-c(testSet2,testSet1[[k]][[1]])
  }
  testSet<-testSet2
  
  SET1<-strsplit(names(referenceSet$x),"_")
  SET2<-c()
  for(k in 1:length(SET1)){
    SET2<-c(SET2,SET1[[k]][[1]])
  }
  AllGenes<-union(SET2,allGenes)
  
  analyticResult <-foreach (curPathway = allPathways) %dopar%{    
    mSigDB_index <- which(DB$geneset.names == curPathway)
    curPathwayGenes <- DB$genesets[[mSigDB_index]]
    
    # R5 class call 
    # source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R") # this is for BigR
    pathwayAnalysis<-myPathwayAnalysis$new()
    pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
    pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
    print(length(intersect(curPathwayGenes,testSet)))
    return(pathwayAnalysis)
    
  }
  savename=paste(Folder,"BIOCARTA_",kk,".Rdata",sep="")
  save(analyticResult,file=savename)
  #AllResults[[kk]]<-analyticResult
}

allPathways <- mSigDB_annotations$objects$C2$REACTOME

# preparing for FET : make total set
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
allGenes <-c()
for (i in 1:length(geneAllSetList)){
  allGenes<-union(allGenes,geneAllSetList[[i]])
}

Folder <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/NP/Ridge2/"

#AllResults<-list()

for(kk in 1:24){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(Folder,"fsDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  
  MAT<-as.matrix(resultsRidge)
  NAME <- rownames(MAT)
  MAT1<-MAT[-1]
  names(MAT1)<-NAME[-1]
  referenceSet<-sort(abs(MAT1), decreasing =T, index.return =T)
  
  
  testSet = names(referenceSet$x[1:500])
  testSet1 = strsplit(testSet,"_")
  testSet2<-c()
  for(k in 1:length(testSet1)){
    testSet2<-c(testSet2,testSet1[[k]][[1]])
  }
  testSet<-testSet2
  
  SET1<-strsplit(names(referenceSet$x),"_")
  SET2<-c()
  for(k in 1:length(SET1)){
    SET2<-c(SET2,SET1[[k]][[1]])
  }
  AllGenes<-union(SET2,allGenes)
  
  analyticResult <-foreach (curPathway = allPathways) %dopar%{    
    mSigDB_index <- which(DB$geneset.names == curPathway)
    curPathwayGenes <- DB$genesets[[mSigDB_index]]
    
    # R5 class call 
    # source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R") # this is for BigR
    pathwayAnalysis<-myPathwayAnalysis$new()
    pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
    pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
    print(length(intersect(curPathwayGenes,testSet)))
    return(pathwayAnalysis)
    
  }
  savename=paste(Folder,"REACTOME_",kk,".Rdata",sep="")
  save(analyticResult,file=savename)
  #AllResults[[kk]]<-analyticResult
}



allPathways <- mSigDB_annotations$objects$C2$KEGG

# preparing for FET : make total set
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
allGenes <-c()
for (i in 1:length(geneAllSetList)){
  allGenes<-union(allGenes,geneAllSetList[[i]])
}

Folder <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/NP/Ridge2/"

#AllResults<-list()

for(kk in 1:24){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(Folder,"fsDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  
  MAT<-as.matrix(resultsRidge)
  NAME <- rownames(MAT)
  MAT1<-MAT[-1]
  names(MAT1)<-NAME[-1]
  referenceSet<-sort(abs(MAT1), decreasing =T, index.return =T)
  
  
  testSet = names(referenceSet$x[1:500])
  testSet1 = strsplit(testSet,"_")
  testSet2<-c()
  for(k in 1:length(testSet1)){
    testSet2<-c(testSet2,testSet1[[k]][[1]])
  }
  testSet<-testSet2
  
  SET1<-strsplit(names(referenceSet$x),"_")
  SET2<-c()
  for(k in 1:length(SET1)){
    SET2<-c(SET2,SET1[[k]][[1]])
  }
  AllGenes<-union(SET2,allGenes)
  
  analyticResult <-foreach (curPathway = allPathways) %dopar%{    
    mSigDB_index <- which(DB$geneset.names == curPathway)
    curPathwayGenes <- DB$genesets[[mSigDB_index]]
    
    # R5 class call 
    # source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R") # this is for BigR
    pathwayAnalysis<-myPathwayAnalysis$new()
    pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
    pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
    print(length(intersect(curPathwayGenes,testSet)))
    return(pathwayAnalysis)
    
  }
  savename=paste(Folder,"KEGG_",kk,".Rdata",sep="")
  save(analyticResult,file=savename)
  #AllResults[[kk]]<-analyticResult
}

