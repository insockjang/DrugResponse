require(predictiveModeling)
require(synapseClient)
require(affy)
require(RColorBrewer)
require(gplots)
require(ggplot2)
require(graphite)
require(RCytoscape)
require(Rgraphviz)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@",mode = "hmac")
source("~/COMPBIO/trunk/users/jang/R5/myPathwayAnalysis.R")
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  

load("~/COMPBIO/trunk/users/jang/pathway_analysis/graphite_pathways.Rdata")

output<-function(resultsScale){
  pred<-c()
  obs<-c()
  pred1<-c()
  obs1<-c()
  for(k in 1:length(resultsScale)){
    pred<-c(pred,resultsScale[[k]]$testPredictions)
    obs<-c(obs,resultsScale[[k]]$testObservations)
    pred1<-c(pred1,resultsScale[[k]]$trainPredictions)
    obs1<-c(obs1,resultsScale[[k]]$trainObservations)
  }
  #   par(mfrow = c(1,2))
  #   plot(obs,pred)
  #   plot(obs1,pred1)
  #   print(c(cor(obs,pred),cor(obs1,pred1)))
  return(cor(obs,pred))
}


cvPredictivePerformance<-function(dirNames,kk,prior = FALSE){
  
  indTestCor<-c()
  indTrainCor<-c()
  overallTestCor<-c()
  overallTrainCor<-c()
  
  for(kk in 1:kk){
    
    filename = paste(dirNames,"/cvDrug_",kk,".Rdata",sep="")
    if(is.element(prior,"CGC")){
      resultsScale <- resultsCGC      
    }
    if(is.element(prior,"TF")){
      resultsScale <- resultsTF      
    }
    if(is.element(prior,"M")){
      resultsScale <- resultsM      
    }
    if(is.element(prior,"R")){
      resultsScale <- resultsR      
    }
    
    load(filename)
    
    testCor <- c()
    trainCor <- c()
    
    for(k in 1:length(resultsScale)){
      testCor <- c(testCor,cor(resultsScale[[k]]$testPredictions,resultsScale[[k]]$testObservations))
      trainCor <- c(testCor,cor(resultsScale[[k]]$trainPredictions,resultsScale[[k]]$trainObservations))
    }
    
    indTestCor<-cbind(indTestCor,testCor)
    indTrainCor<-cbind(indTrainCor,trainCor)
    
    trPred <- foreach(k = 1:length(resultsScale)) %do%{resultsScale[[k]]$trainPredictions}
    tePred <- foreach(k = 1:length(resultsScale)) %do%{resultsScale[[k]]$testPredictions}
    trObsr <- foreach(k = 1:length(resultsScale)) %do%{resultsScale[[k]]$trainObservations}
    teObsr <- foreach(k = 1:length(resultsScale)) %do%{resultsScale[[k]]$testObservations}
    
    allTrPred<-do.call("c",trPred)
    allTePred<-do.call("c",tePred)
    allTrObsr<-do.call("c",trObsr)
    allTeObsr<-do.call("c",teObsr)
    
    allTestCor<-cor(allTePred,allTeObsr)
    allTrainCor<-cor(allTrPred,allTrObsr)
    
    overallTestCor<-c(overallTestCor,allTestCor)
    overallTrainCor<-c(overallTrainCor,allTrainCor)
  }
  return(list(indTestCor=indTestCor,indTrainCor=indTrainCor,overallTestCor=overallTestCor,overallTrainCor=overallTrainCor))
}


pvalMatrix<-function(file,pathwayDBName){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # CCLE drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  # MSigDB from synapse
  mSigDB_annotations <- loadEntity(105363)
  mSigDB_symbolID <- loadEntity(105350)
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID
  if(is.element(pathwayDBName,"KEGG")) allPathways <- mSigDB_annotations$objects$C2$KEGG
  if(is.element(pathwayDBName,"BIOCARTA")) allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  if(is.element(pathwayDBName,"REACTOME")) allPathways <- mSigDB_annotations$objects$C2$REACTOME
  if(is.element(pathwayDBName,"GO_BP")) allPathways <- mSigDB_annotations$objects$C5$GO_BP
  if(is.element(pathwayDBName,"GO_CC")) allPathways <- mSigDB_annotations$objects$C5$GO_CC
  if(is.element(pathwayDBName,"GO_MF")) allPathways <- mSigDB_annotations$objects$C5$GO_MF
  
  pval<-c()
  nes<-c()
  for(k in 1:length(drugName)){
    #   filename = paste("~/PredictiveModel/drugResponse_woCopy/CCLE/ENet2/REACTOME_",k,".Rdata",sep = "")
    filename = paste(file,"/",pathwayDBName,"_",k,".Rdata",sep = "")
    load(filename)
    c1<-c()
    c2<-c()
    for(k1 in 1:length(analyticResult)){
      c1<-c(c1,analyticResult[[k1]]$gseaResult$p.value)
      c2<-c(c2,analyticResult[[k1]]$gseaResult$nes)
    }
    pval<-cbind(pval,c1)
    nes<-cbind(nes,c2)
  }
  colnames(pval)<-drugName
  colnames(nes)<-drugName
  rownames(pval)<-allPathways
  rownames(nes)<-allPathways
  #
  #   c1<-c()
  #   c2<-c()
  #   c3<-c()
  #   for(k in 1:24){
  #     c1<-c(c1,which.min(pval1[,k]))
  #     c2<-c(c2,which.min(pval2[,k]))
  #     c3<-c(c3,which.min(pval3[,k]))
  #   }
  return(list(PVAL = pval, NES = nes))
}


newPvalMatrix<-function(file,pathwayDBName){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # CCLE drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  # MSigDB from synapse
  
  if(is.element(pathwayDBName,"KEGG")){
    allPathways <-c()
    for(k in 1:length(kegg)){
      allPathways <-c(allPathways,kegg[[k]]@title)
    }
  }
  if(is.element(pathwayDBName,"BIOCARTA")){
    allPathways <-c()
    for(k in 1:length(biocarta)){
      allPathways <-c(allPathways,biocarta[[k]]@title)
    }
  }
  if(is.element(pathwayDBName,"REACTOME")){
    allPathways <-c()
    for(k in 1:length(reactome)){
      allPathways <-c(allPathways,reactome[[k]]@title)
    }
  }
  if(is.element(pathwayDBName,"NCI")){
    allPathways <-c()
    for(k in 1:length(nci)){
      allPathways <-c(allPathways,nci[[k]]@title)
    }
  }
  
  
  pval<-c()
  nes<-c()
  for(k in 1:length(drugName)){
    #   filename = paste("~/PredictiveModel/drugResponse_woCopy/CCLE/ENet2/REACTOME_",k,".Rdata",sep = "")
    filename = paste(file,"/",pathwayDBName,"_graphite_",k,".Rdata",sep = "")
    load(filename)
    c1<-c()
    c2<-c()
    for(k1 in 1:length(analyticResult)){
      c1<-c(c1,analyticResult[[k1]]$gseaResult$p.value)
      c2<-c(c2,analyticResult[[k1]]$gseaResult$nes)
    }
    pval<-cbind(pval,c1)
    nes<-cbind(nes,c2)
  }
  colnames(pval)<-drugName
  colnames(nes)<-drugName
  rownames(pval)<-make.names(allPathways)
  rownames(nes)<-make.names(allPathways)
  #
  #   c1<-c()
  #   c2<-c()
  #   c3<-c()
  #   for(k in 1:24){
  #     c1<-c(c1,which.min(pval1[,k]))
  #     c2<-c(c2,which.min(pval2[,k]))
  #     c3<-c(c3,which.min(pval3[,k]))
  #   }
  return(list(PVAL = pval, NES = nes))
}



fetMatrix<-function(file,pathwayDBName){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # CCLE drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  # MSigDB from synapse
  mSigDB_annotations <- loadEntity(105363)
  mSigDB_symbolID <- loadEntity(105350)
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID
  if(is.element(pathwayDBName,"KEGG")) allPathways <- mSigDB_annotations$objects$C2$KEGG
  if(is.element(pathwayDBName,"BIOCARTA")) allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  if(is.element(pathwayDBName,"REACTOME")) allPathways <- mSigDB_annotations$objects$C2$REACTOME
  if(is.element(pathwayDBName,"GO_BP")) allPathways <- mSigDB_annotations$objects$C5$GO_BP
  if(is.element(pathwayDBName,"GO_CC")) allPathways <- mSigDB_annotations$objects$C5$GO_CC
  if(is.element(pathwayDBName,"GO_MF")) allPathways <- mSigDB_annotations$objects$C5$GO_MF
  
  pval<-c()
  oddRatio<-c()
  for(k in 1:length(drugName)){
    
    filename = paste(file,"/",pathwayDBName,"_",k,".Rdata",sep = "")
    load(filename)
    c1<-c()
    c2<-c()
    for(k1 in 1:length(analyticResult)){
      c1<-c(c1,analyticResult[[k1]]$fetResult$p.value)
      c2<-c(c2,analyticResult[[k1]]$fetResult$estimate)
    }
    pval<-cbind(pval,c1)
    oddRatio<-cbind(oddRatio,c2)
  }
  colnames(pval)<-drugName
  colnames(oddRatio)<-drugName
  rownames(pval)<-allPathways
  rownames(oddRatio)<-allPathways
  
  return(list(PVAL = pval, ODDRATIO = oddRatio))
}


AUCMatrix<-function(file,method = "median_mad"){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # Sanger drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  auc<-c()
  for(k in 1:length(drugName)){
    filename = paste(file,"/ROCR_",method,"_",k,".Rdata",sep = "")
    load(filename)    
    auc<-c(auc,resultsCat$erAUC@y.values[[1]])  
  }
  names(auc)<-drugName
  
  return(list(AUC = auc))
}


corMatrix<-function(file){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # CCLE drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  performance<-c()
  for(k in 1:length(drugName)){
    filename = paste(file,"/cvDrug_",k,".Rdata",sep = "")
    load(filename)
    performance<-c(performance,output(resultsScale))    
  }
  names(performance)<-drugName
  
  return(list(COR = performance))
}



pvalSigMatrix<-function(file,pathwayDBName,geneNum = 0){
  if(length(grep("CCLE",file))==1){
    id_drugLayer <- "269024" # CCLE drugs 24
  }
  if(length(grep("Sanger",file))==1){
    id_drugLayer <- "220680" # CCLE drugs 24
  }
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  # MSigDB from synapse
  mSigDB_annotations <- loadEntity(105363)
  mSigDB_symbolID <- loadEntity(105350)
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID
  if(is.element(pathwayDBName,"KEGG")) allPathways <- mSigDB_annotations$objects$C2$KEGG
  if(is.element(pathwayDBName,"BIOCARTA")) allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  if(is.element(pathwayDBName,"REACTOME")) allPathways <- mSigDB_annotations$objects$C2$REACTOME
  if(is.element(pathwayDBName,"GO_BP")) allPathways <- mSigDB_annotations$objects$C5$GO_BP
  if(is.element(pathwayDBName,"GO_CC")) allPathways <- mSigDB_annotations$objects$C5$GO_CC
  if(is.element(pathwayDBName,"GO_MF")) allPathways <- mSigDB_annotations$objects$C5$GO_MF
  pval<-c()
  nes<-c()
  hit<-c()
  for(k in 1:length(drugName)){
    
    filename = paste(file,"/",pathwayDBName,"_",k,".Rdata",sep = "")
    load(filename)
    c1<-c()
    c2<-c()
    S<-c()
    K<-c()
    pval1<-matrix(1,nrow = length(analyticResult),ncol =1)
    nes1<-matrix(NaN,nrow = length(analyticResult),ncol =1)
    for(k1 in 1:length(analyticResult)){
      c1<-c(c1,analyticResult[[k1]]$gseaResult$p.value)
      c2<-c(c2,analyticResult[[k1]]$gseaResult$nes)
      S<-c(S,sum(analyticResult[[k1]]$gseaResult$hit))      
    }
    a1 <- which(c1<=0.05)
    a2 <- which(S>=geneNum)
    a3 <- intersect(a1,a2)
    pval1[a3]<-c1[a3]
    nes1[a3]<-c2[a3]
    
    pval<-cbind(pval,pval1)
    nes<-cbind(nes,nes1)
    hit<-cbind(hit,S)
  }
  colnames(pval)<-drugName
  colnames(nes)<-drugName
  colnames(hit)<-drugName
  rownames(pval)<-allPathways
  rownames(nes)<-allPathways
  rownames(hit)<-allPathways
  #
  #   c1<-c()
  #   c2<-c()
  #   c3<-c()
  #   for(k in 1:24){
  #     c1<-c(c1,which.min(pval1[,k]))
  #     c2<-c(c2,which.min(pval2[,k]))
  #     c3<-c(c3,which.min(pval3[,k]))
  #   }
  return(list(PVAL = pval, NES = nes, HIT = hit))
}



heatmapPathway<-function(folderName,pathwayName,NES_threshold=0,PVAL_threshold=0.05,HIT_threshold = 0){
  require(ggplot2)
  
  B1<-pvalSigMatrix(folderName,pathwayName)
  
  coeff_nes<- B1$NES > NES_threshold
  coeff_nes[is.na(coeff_nes)]<-0
  coeff_hit<-B1$HIT>= HIT_threshold
  coeff_hit[is.na(coeff_hit)]<-0
  
  coeff_pval<-B1$PVAL<= PVAL_threshold
  
  coeff <- coeff_nes * coeff_hit * coeff_pval
  b2<-B1$PVAL
  
  b2[is.na(b2)]<-1
  pp<-b2
  pp[which(pp==0)]<-10^-4
  
  processed<--log10(pp)*coeff
  
  aa<-apply(processed,1,sum)
  aaa<-which(aa!=0)
  bb<-apply(processed,2,sum)
  bbb<-which(bb!=0)
  
  allPathways<-rownames(B1$NES)
  a<-heatmap.2(processed[aaa,bbb], Rowv=T, Colv=T, 
               col = colorRampPalette(brewer.pal(5,"BrBG"))(256), scale="none", key=T, density.info="none", 
               trace="none", cexRow=1, cexCol=1, symm=F,symkey=F,symbreaks=T)
  return(list(heatmap = a, pathway = allPathways[aaa]))
  
}


pathwaySynapse<-function(file,pathwayDB,parentID){
  ENet<-pvalMatrix(paste(file,"/ENet2",sep=""),pathwayDB)
  Lasso<-pvalMatrix(paste(file,"/Lasso2",sep=""),pathwayDB)
  Ridge<-pvalMatrix(paste(file,"/Ridge2",sep=""),pathwayDB)
  bsENet<-pvalMatrix(paste(file,"/bootstrapENet2",sep=""),pathwayDB)
  bsLasso<-pvalMatrix(paste(file,"/bootstrapLasso2",sep=""),pathwayDB)
  bsRidge<-pvalMatrix(paste(file,"/bootstrapRidge2",sep=""),pathwayDB)
  
  
  ## REQUIRED VALUES ARE name, parentId, AND type
  myData1 <- ExpressionData(list(name=paste(pathwayDB,"_ENet",sep=""),parentId = parentID))
  myData1 <- createEntity(myData1)
  myData1 <- addObject(myData1, ENet)
  storeEntity(myData1)
  
  myData2 <- ExpressionData(list(name=paste(pathwayDB,"_Lasso",sep=""),parentId = parentID))
  myData2 <- createEntity(myData2)
  myData2 <- addObject(myData2, Lasso)
  storeEntity(myData2)
  
  myData3 <- ExpressionData(list(name=paste(pathwayDB,"_Ridge",sep=""),parentId = parentID))
  myData3 <- createEntity(myData3)
  myData3 <- addObject(myData3, Ridge)
  storeEntity(myData3)
  
  myData4 <- ExpressionData(list(name=paste(pathwayDB,"_boostrapENet",sep=""),parentId = parentID))
  myData4 <- createEntity(myData4)
  myData4 <- addObject(myData4, bsENet)
  storeEntity(myData4)
  
  myData5 <- ExpressionData(list(name=paste(pathwayDB,"_bootstrapLasso",sep=""),parentId = parentID))
  myData5 <- createEntity(myData5)
  myData5 <- addObject(myData5, bsLasso)
  storeEntity(myData5)
  
  myData6 <- ExpressionData(list(name=paste(pathwayDB,"_bootstrapRidge",sep=""),parentId = parentID))
  myData6 <- createEntity(myData6)
  myData6 <- addObject(myData6, bsRidge)
  storeEntity(myData6)
}


hclusterAnalysis<-function(pENet,Title){
  
  pENet[is.na(pENet)]<-1
  
  # coeff<- pENet<=0.05
  
  pp<-pENet
  pp[which(pp==0)]<-10^-4
  
  # processed<--log10(pp)*coeff
  processed<--log10(pp)
  
  # png("HC_all_bsENet.png",width=800,height=800)
  a<-heatmap.2(processed, Rowv=T, Colv=T, 
               col = colorRampPalette(brewer.pal(5,"BrBG"))(256), scale="none", key=T, density.info="none", 
               trace="none", cexRow=1, cexCol=1, symm=F,symkey=F,symbreaks=T,main = Title)
  return(a)
}

fetSynapse<-function(file,pathwayDB,parentID){
  ENet<-fetMatrix(paste(file,"/ENet2",sep=""),pathwayDB)
  Lasso<-fetMatrix(paste(file,"/Lasso2",sep=""),pathwayDB)
  Ridge<-fetMatrix(paste(file,"/Ridge2",sep=""),pathwayDB)
  bsENet<-fetMatrix(paste(file,"/bootstrapENet2",sep=""),pathwayDB)
  bsLasso<-fetMatrix(paste(file,"/bootstrapLasso2",sep=""),pathwayDB)
  bsRidge<-fetMatrix(paste(file,"/bootstrapRidge2",sep=""),pathwayDB)
  
  
  ## REQUIRED VALUES ARE name, parentId, AND type
  myData1 <- ExpressionData(list(name=paste(pathwayDB,"_ENet",sep=""),parentId = parentID))
  myData1 <- createEntity(myData1)
  myData1 <- addObject(myData1, ENet)
  storeEntity(myData1)
  
  myData2 <- ExpressionData(list(name=paste(pathwayDB,"_Lasso",sep=""),parentId = parentID))
  myData2 <- createEntity(myData2)
  myData2 <- addObject(myData2, Lasso)
  storeEntity(myData2)
  
  myData3 <- ExpressionData(list(name=paste(pathwayDB,"_Ridge",sep=""),parentId = parentID))
  myData3 <- createEntity(myData3)
  myData3 <- addObject(myData3, Ridge)
  storeEntity(myData3)
  
  myData4 <- ExpressionData(list(name=paste(pathwayDB,"_boostrapENet",sep=""),parentId = parentID))
  myData4 <- createEntity(myData4)
  myData4 <- addObject(myData4, bsENet)
  storeEntity(myData4)
  
  myData5 <- ExpressionData(list(name=paste(pathwayDB,"_bootstrapLasso",sep=""),parentId = parentID))
  myData5 <- createEntity(myData5)
  myData5 <- addObject(myData5, bsLasso)
  storeEntity(myData5)
  
  myData6 <- ExpressionData(list(name=paste(pathwayDB,"_bootstrapRidge",sep=""),parentId = parentID))
  myData6 <- createEntity(myData6)
  myData6 <- addObject(myData6, bsRidge)
  storeEntity(myData6)
}


pathwayFet<-function(allPathways){
  
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
  allGenes <-c()
  for (i in 1:length(geneAllSetList)){
    allGenes<-union(allGenes,geneAllSetList[[i]])
  }
  
  pVal<-matrix(0,nrow = length(allPathways),ncol = length(allPathways))
  oddRatio<-matrix(0,nrow = length(allPathways),ncol = length(allPathways))
  
  for(k1 in 1:length(allPathways)){
    geneSet1<-DB$genesets[is.element(DB$geneset.names,allPathways[k1])][[1]]
    for(k2 in k1:length(allPathways)){
      geneSet2<-DB$genesets[is.element(DB$geneset.names,allPathways[k2])][[1]]
      
      Mat2x2 <- mat.or.vec(2,2)
      Mat2x2[1,1] <- length(intersect(geneSet1,geneSet2))
      Mat2x2[2,1] <- length(setdiff(geneSet1,geneSet2))
      Mat2x2[1,2] <- length(setdiff(geneSet2,geneSet1))
      Mat2x2[2,2] <- length(allGenes) - Mat2x2[1,1]- Mat2x2[1,2]- Mat2x2[2,1]
      
      fet<-fisher.test(Mat2x2)                            
      pVal[k1,k2]<-fet$p.value
      oddRatio[k1,k2]<-fet$estimate
    }
  }
  pValU<-triu(pVal,1)
  pValL<-t(pVal)
  pVal<-pValU+pValL+diag(diag(pVal))
  rownames(pVal)<-allPathways
  colnames(pVal)<-allPathways
  
  return(pVal)
  
}


convertPathwayGraph<-function(pathwayObject,isolatedNode){
  
  Node<-nodes(pathwayObject)
  Edge<-edges(pathwayObject)
  
  # initialize
  g = new ("graphNEL", edgemode = "directed")   
  g = initNodeAttribute (g, "nodeType", "char", "undefined")
  g = initNodeAttribute (g, "label", "char", "undefined")
  g = initEdgeAttribute (g, "edgeType", "char", "undefined")
  
  if(isolatedNode){
    #   add nodes
    for(k in 1:length(Node)){
      g = addNode(Node[k],g)
      nodeData (g, Node[k], 'nodeType') = Node[k]
    }
    
    if(nrow(Edge)>0){
      edge_processed<-unique(Edge[,c("src","dest")])    
      
      
      #   add edges
      for(k in 1:nrow(edge_processed)){
        g = addEdge(edge_processed$src[k],edge_processed$dest[k],g)
        edgeData(g, edge_processed[k,1],edge_processed[k,2],"edgeType") = "directed"
      }
    }
  }
  else{
    if(nrow(Edge)>0){
      edge_processed<-unique(Edge[,c("src","dest")])    
      
      Node1<-union(edge_processed[,1],edge_processed[,2])
      
      #   add nodes
      for(k in 1:length(Node1)){
        g = addNode(Node1[k],g)
        nodeData (g, Node1[k], 'nodeType') = Node1[k]
      }
      #   add edges
      for(k in 1:nrow(edge_processed)){
        g = addEdge(edge_processed$src[k],edge_processed$dest[k],g)
        edgeData(g, edge_processed[k,1],edge_processed[k,2],"edgeType") = "directed"
      }
    }
  }
  return(g)
}

# 
# pathwayPlot<-function(graphObject,features,outType="dot",filename= "temp"){
#   
#   g<-graphObject    
#   names(labels) = labels = nodes(g)
#   namePathway <- names(labels)
#   
#   if(!is.nan(features)){     
#     # selected feature 
#     Name<-names(features)
#     X<-sort(features,index.return=TRUE)
#     Name<-names(X$x)
#     
#     Name1<-strsplit(Name,"_")  
#     Name2<-c()
#     for(k in 1:length(Name)){
#       Name2<-c(Name2,Name1[[k]][1])
#     }  
#     
#     Colors<-colorRampPalette(brewer.pal(9,"Greens"))(length(features))
#     
#     b<-match(namePathway,Name2)
#     b[is.na(b)]<-0
#   }
#   else{
#     b<-sapply(labels, function(x) {0})
#   }
#   fontscale = round(200/sqrt(length(labels)))
# #   fontscale = round(1400/length(labels))
#   #   fontsize = fontscale * max(nchar(unlist(strsplit(labels,"\n", fixed=TRUE))))
#   #   fontsize = fontscale
#   
#   fill <- sapply(labels, function(x) {0})
#   for(k in 1:length(fill)){
#     if(b[k]==0){
#       fill[k] = "yellow"
#     }
#     else{
#       fill[k]=Colors[b[k]]
#     }
#   }
#   shape <- sapply(labels, function(x) {"circle"})  
#   width<-sapply(labels, function(x) {10})
#   height<-sapply(labels, function(x) {10})
#   fontsize<-sapply(labels, function(x) {20})
#   
#   names(fill) = names(shape) = names(fontsize) = nodes(g)
#   
#   # Note: fill, shape, width & height are all vectors which I set to get custom
#   # node shapes, sizes, and fill color.
#   attrs <- list(graph=list(rankdir="LR"))
#   attrs <- list(edge = list(arrowsize = 1.0))
#   attrs <- list(node = list(fixedsize = TRUE))
#   
#   x <- layoutGraph(g,layoutType = outType,attrs = attrs) # “dot”, “neato”, “twopi”, “circo”, “fdp”
#   
#   nodeRenderInfo(x) = list(label = labels, shape = shape, fontsize = fontsize, fill = fill,height = height,width = width)
#   
# #   png(paste(filename,".png",sep=""),width = 600*(round(length(labels)/20)+1),height=400*(round(length(labels)/20)+1))
#   renderGraph(x)  
#   #   title(main = )
# #   dev.off()
#   #   return(x)
#   
# }
