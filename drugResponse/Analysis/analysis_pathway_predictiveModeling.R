library(affy)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
# MSigDB from synapse
mSigDB_annotations <- loadEntity(105363)
mSigDB_symbolID <- loadEntity(105350)
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  
allPathways <- mSigDB_annotations$objects$C2$KEGG


id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF

drugName<-colnames(pData(adf_drug))

##################################################################################################################
# upload synapse 
##################################################################################################################
myName <- "In Sock Jang"
project <- Project(list(
  name=paste("CCLE predictiveModeling - ", myName)
  ))
project <- createEntity(project)

onWeb(project)

dataset <- Dataset(list(
  name="Lasso Bootstrap Pathway Analysis",
  parentId="161832"
  ))
(dataset <- createEntity(dataset))
datasetID<-propertyValue(dataset, "id")
##################################################################################################################
##################################################################################################################
ids<-c()

for(k in 1:24){
  filename = paste("Lasso_original/lassoFeatureAnalysisWithBootstrap_",k,".Rdata",sep ="")
  load(filename)
  
  ref<-referenceSet[-1]
  names(ref) <- rownames(referenceSet)[-1]
  
  a<-sort(abs(ref), decreasing = TRUE, index.return = TRUE)
  
  referenceGenes<-names(ref)[a$ix[1:200]]
  getUniqueGenesFromFeatureNames <- function(featureNames){
    genes_copy <- sub("_copy","",featureNames[grep("_copy", featureNames)])
    genes_expr <- sub("_expr","",featureNames[grep("_expr", featureNames)])
    genes_mut <- sub("_mut","",featureNames[grep("_mut", featureNames)])  
    geneSet <- union(genes_copy,union(genes_expr, genes_mut))
  }
  
  referenceGenes<-getUniqueGenesFromFeatureNames(referenceGenes)
  
  gsea.pval<-c()
  for(kk in 1:length(analysis)){
    if(is.null(analysis[[kk]]$GSEA$p.value)){
      gsea.pval<-c(gsea.pval,1)
    }
    else{
      gsea.pval<-c(gsea.pval,analysis[[kk]]$GSEA$p.value)    
    }        
  }
  filename1 = paste("Lasso_original/lassoFETAnalysisWithBootstrap_",k,".Rdata",sep ="")
  load(filename1)
  
  fet.pval<-c()
  for(kk in 1:length(FET)){
    if(is.null(FET[[k]]$p.value)){
      fet.pval<-c(fet.pval,1)
    }
    else{
      fet.pval<-c(fet.pval,FET[[kk]]$p.value)    
    }
  }
  
  statResults =list(GSEA = gsea.pval,FET = fet.pval)
  
  objnames <- paste("DRUG_",drugName[k],"_pathwayAnalysis")
  
  submittedModelDatasetId <- datasetID
  submittedModelLayer <- Layer(list(name = objnames, type = "E", parentId = submittedModelDatasetId))
  submittedModelLayer <- addObject(entity=submittedModelLayer, object=statResults)
  submittedModelLayer <- storeEntity(submittedModelLayer)
  ids<-c(ids,propertyValue(submittedModelLayer, "id"))
}
#synapseIDs_drug <- matrix(0,nrow = 24,ncol = 5)
synapseIDs_drug[,1]=ids
