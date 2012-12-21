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

refList<-as.numeric(mat.or.vec(length(set),1))
names(refList)<-set

caretModel_eNet <- predictiveModel_eNet$rawCaretModel()
coefs_eNet <- caretModel_eNet$finalModel$beta[, ncol(caretModel_eNet$finalModel$beta)]
featureNames_sig <- names(coefs_eNet[coefs_eNet != 0] )
feature_sig <- coefs_eNet[coefs_eNet != 0]

for(i in 1:length(feature_sig)){
  if(length(grep("_expr",names(feature_sig)[i])))
    names(feature_sig[i]) <- sub("_expr", "", featureNames_sig[i])
  if(length(grep("_copy",names(feature_sig)[i])))
    names(feature_sig[i]) <- sub("_copy", "", featureNames_sig[i])
  if(length(grep("_mut",names(feature_sig)[i])))
    names(feature_sig[i]) <- sub("_mut", "", featureNames_sig[i])
}

if (length(sig_set) != length(sig_set1) + length(sig_set1) + length(sig_set1)){
  print("ERROR")
  print("Figure out which feature should be selected among the same features with different coefficients")
}
####################################################
#### In reference list, I included all feature from featuredata
###################################################
# assign the coefficient to newly generated reference list 
refList[names(refList) %in% names(feature_sig)]<-feature_sig[match(names(refList),names(feature_sig))[which(!is.na(match(names(refList),names(feature_sig))))]]

mSigDB_annotations <- loadEntity(105363)
mSigDB_symbolID <- loadEntity(105350)

DB<-mSigDB_symbolID$objects$MsigDB_symbolID
allPathways <- mSigDB_annotations$objects$C2$BIOCARTA

pval <-c()
es<-c()


for (curPathway in allPathways){
  print(paste("processing", curPathway))
  
  mSigDB_index <- which(DB$geneset.names == curPathway)
  curPathwayGenes <- DB$genesets[[mSigDB_index]]
  
  #if (length(setdiff(curPathwayGenes,names(refList)))){
  #  print("error")
  #  break
  #}
  curPathwayGenes <- intersect(curPathwayGenes, names(refList))
  
  gseaResults <- preRankedTest(abs(refList), curPathwayGenes)
  pval <-rbind(pval, gseaResults$Pval)
  es <- rbind(es,gseaResults$ES)
}

getUniqueGenesFromFeatureNames <- function(featureNames){
  genes_copy <- featureNames[grep("_copy", featureNames)]
  genes_expr <- featureNames[grep("_expr", featureNames)]
  genes_mut <- featureNames[grep("_mut", featureNames)]
  
  geneSet <- union(genes_copy,union(genes_expr, genes_mut))
}