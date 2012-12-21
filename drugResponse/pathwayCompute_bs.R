library(synapseClient)
library(predictiveModeling)
library(multicore)
library(doMC)
registerDoMC()

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  


bootstrapFolder <- "~/COMPBIO/trunk/users/jang/drugResponse/Sanger/NP/bootstrapRidge2/"

allPathways <- mSigDB_annotations$objects$C2$REACTOME

# preparing for FET : make total set
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
allGenes <-c()
for (i in 1:length(geneAllSetList)){
  allGenes<-union(allGenes,geneAllSetList[[i]])
}

for(kk in 1:130){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(bootstrapFolder,"cvDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  resultsBsModel<-resultsRidge
  
  # processing to find significant bootstrapping feature selection
  ResultBS<-c()
  for(k in 1:length(resultsBsModel)){    
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsBsModel[[k]][-1])),ties.method="min")/length(resultsBsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsBsModel[[k]])[-1]  
  reference <- apply(ResultBS,1,sum)
  
  referenceSet<-sort(reference, decreasing =T, index.return =T)
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
  savename=paste(bootstrapFolder,"REACTOME_",kk,".Rdata",sep="")
  save(analyticResult,file=savename)
  #AllResults[[kk]]<-analyticResult
}


allPathways <- mSigDB_annotations$objects$C2$BIOCARTA

# preparing for FET : make total set
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
allGenes <-c()
for (i in 1:length(geneAllSetList)){
  allGenes<-union(allGenes,geneAllSetList[[i]])
}

for(kk in 1:130){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(bootstrapFolder,"cvDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  resultsBsModel<-resultsRidge
  
  # processing to find significant bootstrapping feature selection
  ResultBS<-c()
  for(k in 1:length(resultsBsModel)){    
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsBsModel[[k]][-1])),ties.method="min")/length(resultsBsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsBsModel[[k]])[-1]  
  reference <- apply(ResultBS,1,sum)
  
  referenceSet<-sort(reference, decreasing =T, index.return =T)
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
  savename=paste(bootstrapFolder,"BIOCARTA_",kk,".Rdata",sep="")
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

for(kk in 1:130){
  
  # this part can be replaced with R5 class for each bootstrapping model
  fileLasso = paste(bootstrapFolder,"cvDrug_",kk,".Rdata",sep="")  
  load(fileLasso)
  resultsBsModel<-resultsRidge
  
  # processing to find significant bootstrapping feature selection
  ResultBS<-c()
  for(k in 1:length(resultsBsModel)){    
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsBsModel[[k]][-1])),ties.method="min")/length(resultsBsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsBsModel[[k]])[-1]  
  reference <- apply(ResultBS,1,sum)
  
  referenceSet<-sort(reference, decreasing =T, index.return =T)
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
  savename=paste(bootstrapFolder,"KEGG_",kk,".Rdata",sep="")
  save(analyticResult,file=savename)
  #AllResults[[kk]]<-analyticResult
}



