MsigDBPerformanceTest <- function(coefs_eNet,MsigDB_SET,DB){
  ### how to use MsigDB for geneset Fisher Exact Test
  ## Fist, specify your geneset
  #geneSet <- "your geneset"
  
  NAME <-names(coefs_eNet)
  NAME1<-NAME[coefs_eNet !=0]
  
  COPY <- NAME1[grep("_copy",NAME1)]
  MUT <- NAME1[grep("_mut",NAME1)]
  EXP <- NAME1[grep("_exp",NAME1)]
  
  COPY1 <- sub("_copy","",COPY)
  MUT1 <- sub("_mut","",MUT)
  EXP1 <- sub("_exp","",EXP)
  
  geneSet <- union(COPY1,union(EXP1,MUT1))
  
  ## Should figure out how to import MsigDB from Synapse
  # DB <- loadEntity() # MsigDB
  # DB$geneset.names
  
  # Annotation <- loadEntity() # MsigDB annotation
  # Annotation$C1$POS
  
  ## Second, specify MsigDB. 
  # e.g. cannonical pathway(KEGG, BiOCARTA, or REACTOME) or Gene Ontology(GO) biological process, cellular component, or molecular function
  
  # Chromosomal Position 
  #   C1$POS
  # Curated Gene Sets
  #   Chemical and genetic perturbation
  #     C2$CGP
  #   Cannonical Pathway: BIOCARTA/KEGG/REACTOME
  #     C2$BIOCARTA
  #     C2$KEGG
  #     C2$REACTOME
  #   Motif : microRNA target / transcription factor target
  #     C3$miRT
  #     C3$TFT
  #   Computational gene set : cancer gene module / cancer gene neighbor
  #     C4$CGM
  #     C4$CGN
  #   Gene Onltology(GO)
  #     C5$GO_BP
  #     C5$GO_CC
  #     C5$GO_MF
  # MsigDB_SET <- C5$GO_BP # You can change this term to "BIOCARTA", "REACTOME", etc.
  
  # find how many sets as overall set ()
  numSet <- length(MsigDB_SET)
  
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,MsigDB_SET)]
  
  geneAllSet <-c()
  
  for (i in 1:numSet){
    geneAllSet<-union(geneAllSet,geneAllSetList[[i]])
  }
  
  pval <-c()
  # Fisher Exact Test
  for (i in 1:numSet){
    Mat2x2 <- mat.or.vec(2,2)
    Mat2x2[1,1] <- length(intersect(geneSet,geneAllSetList[[i]]))
    Mat2x2[2,1] <- length(setdiff(geneSet,geneAllSetList[[i]]))
    Mat2x2[1,2] <- length(setdiff(geneAllSetList[[i]],geneSet))
    Mat2x2[2,2] <- length(union(geneAllSet,geneSet)) - Mat2x2[1,1]- Mat2x2[1,2]- Mat2x2[2,1]
    testResult<-fisher.test(Mat2x2)
    pval<-rbind(pval,testResult$p.value)
  }
  PVAL <- sort(pval,index.return = TRUE)
  pval1 <-PVAL$x 
  names(pval1)<-MsigDB_SET[PVAL$ix]
  
  return(pval1)
}
