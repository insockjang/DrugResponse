source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R")
# pathwayAnalysisFuctions for bootstrapping and nonbootstrapping
pathwayAnalysisBootstrap<-function(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE){
  
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- mSigDB_annotations$objects$C2$KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- mSigDB_annotations$objects$C2$REACTOME
  }
  
  
  # preparing for FET : make total set
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
  allGenes <-c()
  for (i in 1:length(geneAllSetList)){
    allGenes<-union(allGenes,geneAllSetList[[i]])
  }
  
  #   bootstrapFolder <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/NP/bootstrapRidge2/"
  
  for(kk in 1:numCompounds){
    
    # this part can be replaced with R5 class for each bootstrapping model
    fileInput = paste(bootstrapFolder,"cvDrug_",kk,".Rdata",sep="")  
    load(fileInput)
    if(!prior){
      if(is.element(modelName, "Ridge")){
        resultsBsModel<-resultsRidge        
      }
      if(is.element(modelName, "Lasso")){
        resultsBsModel<-resultsLasso
      }
      if(is.element(modelName, "ENet")){
        resultsBsModel<-resultsENet
      }         
      savename=paste(bootstrapFolder,pathwayName,"_",kk,".Rdata",sep="")
    }
    else{
      if(is.element(modelName, "CGC")){
        resultsBsModel<-resultsCGC        
      }
      if(is.element(modelName, "TF")){
        resultsBsModel<-resultsTF
      }
      if(is.element(modelName, "M")){
        resultsBsModel<-resultsM
      }
      if(is.element(modelName, "R")){
        resultsBsModel<-resultsR
      }
      savename=paste(bootstrapFolder,pathwayName,"_",modelName,"_",kk,".Rdata",sep="")
    }
    
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
      
      pathwayAnalysis<-myPathwayAnalysis$new()
      pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
      pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
      #print(length(intersect(curPathwayGenes,testSet)))
      return(pathwayAnalysis)      
    }        
    save(analyticResult,file=savename)
    
  }
}

pathwayAnalysisNonbootstrap<-function(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE){
  
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- mSigDB_annotations$objects$C2$KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- mSigDB_annotations$objects$C2$REACTOME
  }
  
  
  # preparing for FET : make total set
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
  allGenes <-c()
  for (i in 1:length(geneAllSetList)){
    allGenes<-union(allGenes,geneAllSetList[[i]])
  }
  
  #   bootstrapFolder <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/NP/bootstrapRidge2/"
  
  for(kk in 1:numCompounds){
    
    # this part can be replaced with R5 class for each bootstrapping model
    fileInput = paste(bootstrapFolder,"fsDrug_",kk,".Rdata",sep="")  
    load(fileInput)
    if(!prior){
      if(is.element(modelName, "Ridge")){
        resultsBsModel<-resultsRidge        
      }
      if(is.element(modelName, "Lasso")){
        resultsBsModel<-resultsLasso
      }
      if(is.element(modelName, "ENet")){
        resultsBsModel<-resultsENet
      }         
      savename=paste(bootstrapFolder,pathwayName,"_",kk,".Rdata",sep="")
    }
    else{
      if(is.element(modelName, "CGC")){
        resultsBsModel<-resultsCGC        
      }
      if(is.element(modelName, "TF")){
        resultsBsModel<-resultsTF
      }
      if(is.element(modelName, "M")){
        resultsBsModel<-resultsM
      }
      if(is.element(modelName, "R")){
        resultsBsModel<-resultsR
      }
      savename=paste(bootstrapFolder,pathwayName,"_",modelName,"_",kk,".Rdata",sep="")
    }
    
    
    MAT<-as.matrix(resultsBsModel)
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
      #print(length(intersect(curPathwayGenes,testSet)))
      return(pathwayAnalysis)
      
    }     
    save(analyticResult,file=savename)
    
  }
}
