library(predictiveModeling)
library(synapseClient)

data(demoData)

ds_features_cn_mut_ccle <- createAggregateFeatureDataSet(list(copy = copySet, mut = oncomapSet))
dataSets_ccleFeaturesCnMut_sangerChems <- createFeatureAndResponseDataList(ds_features_cn_mut_ccle, sangerADF)

featureData_scaled <- t(scale(t(dataSets_ccleFeaturesCnMut_sangerChems$featureData)))
responseData_scaled <- t(scale(t(dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE])))

predictiveModel_eNet <- CaretModel$new(modelType = "glmnet")
predictiveModel_eNet$train(t(featureData_scaled), t(responseData_scaled), tuneGrid=createENetTuneGrid(alphas=1))

featureNames_all <- rownames(featureData_scaled)
set1 <- sub("_expr", "", featureNames_all[grep("_expr",featureNames_all)])
set2 <- sub("_copy", "", featureNames_all[grep("_copy",featureNames_all)])
set3 <- sub("_mut", "", featureNames_all[grep("_mut",featureNames_all)])
set <- union(set1,union(set2,set3))

caretModel_eNet <- predictiveModel_eNet$rawCaretModel()
coefs_eNet <- caretModel_eNet$finalModel$beta[, ncol(caretModel_eNet$finalModel$beta)]
featureNames_sig <- names(coefs_eNet[coefs_eNet != 0] )
feature_sig <- coefs_eNet[coefs_eNet != 0]


getUniqueGenesFromFeatureNames <- function(featureNames){
  genes_copy <- sub("_copy","",featureNames[grep("_copy", featureNames)])
  genes_expr <- sub("_exprs","",featureNames[grep("_expr", featureNames)])
  genes_mut <- sub("_mut","",featureNames[grep("_mut", featureNames)])  
  geneSet <- union(genes_copy,union(genes_expr, genes_mut))
}

feature_sigName <-getUniqueGenesFromFeatureNames(names(feature_sig))


mSigDB_annotations <- loadEntity(105363)
mSigDB_symbolID <- loadEntity(105350)

DB<-mSigDB_symbolID$objects$MsigDB_symbolID

allPathways <- mSigDB_annotations$objects$C2$KEGG

# process for FET 
geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
geneAllSet <-c()
for (i in 1:length(geneAllSetList)){
  geneAllSet<-union(geneAllSet,geneAllSetList[[i]])
}

pval1 <-c()
pval.fet <-c()
nes1<-c()

i <-0
for (curPathway in allPathways){
  i<- i+1
  print(paste("processing", curPathway))
  
  mSigDB_index <- which(DB$geneset.names == curPathway)
  curPathwayGenes <- DB$genesets[[mSigDB_index]]
  
  refList <- mat.or.vec(length(union(curPathwayGenes,feature_sigName)),1)
  names(refList) <- union(curPathwayGenes,feature_sigName)
  refList[names(refList) %in% feature_sigName]<-feature_sig[match(names(refList),feature_sigName)[which(!is.na(match(names(refList),feature_sigName)))]]
  
  #curPathwayGenes <- intersect(curPathwayGenes, names(refList))
  
  gseaResults <- preRankedTest(abs(refList), curPathwayGenes,np = 1000)
  pval1 <-rbind(pval1, gseaResults$Pval)
  nes1 <- rbind(nes1,gseaResults$NES)
  
  Mat2x2 <- mat.or.vec(2,2)
  Mat2x2[1,1] <- length(intersect(feature_sigName,geneAllSetList[[i]]))
  Mat2x2[2,1] <- length(setdiff(feature_sigName,geneAllSetList[[i]]))
  Mat2x2[1,2] <- length(setdiff(geneAllSetList[[i]],feature_sigName))
  Mat2x2[2,2] <- length(union(geneAllSet,feature_sigName)) - Mat2x2[1,1]- Mat2x2[1,2]- Mat2x2[2,1]
  testResult<-fisher.test(Mat2x2)
  pval.fet<-rbind(pval.fet,testResult$p.value)
  
}

allPathways[intersect(which(pval1<=0.05),which(nes1!=0))]
A.gsea<-cbind(t(t(allPathways[which(pval1<=0.05)] )),pval1[which(pval1<=0.05)])
A.fet<-cbind(t(t(allPathways[which(pval.fet<=0.05)])),pval.fet[which(pval.fet <=0.05)])
A.intersect<-intersect(allPathways[which(pval1<=0.05)],allPathways[which(pval.fet<=0.05)])

write.table(A.gsea,file = "gsea_demo.txt",row.name = FALSE,col.name = FALSE)
write.table(A.fet,file = "fet_demo.txt",row.name = FALSE,col.name = FALSE)
write.table(A.intersect,file = "overlap_demo.txt",row.name = FALSE,col.name = FALSE)
