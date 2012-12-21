library(predictiveModeling)
library(affy)
library(synapseClient)
library(RColorBrewer)
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
# synapseID  
##################################################################################################################
# Lasso_bootstrap, Lasso, Enet, Ridge, Spls

ids<-loadEntity("162185")
synapseIDs_drug<-ids$objects$synapseIDs_drug

analysis <- Dataset(list(
  name="New Pathway Analysis Comparison",
  description="GSEA and FET with CCLE Drug Response depending upon predictiveModel",
  parentId="161936"
  ))
(dataset <- createEntity(analysis))

datasetID<-propertyValue(dataset, "id")
##################################################################################################################
# Start analysis comparisons
##################################################################################################################
modelname<-c("Lasso_bootstrap", "Lasso", "Enet", "Ridge", "Spls")
hcResult<-list()
allDrugStatistics<-list()

for(k in 1:5){ # Lasso_boots


  ## select drug
  drugGsea<-c()
  drugFet<-c()
  
  for(kk in 1:24){
    
    print(paste("Drug Name is ", drugName[kk],sep = ""))
    ##
    pathways<-loadEntity(synapseIDs_drug[kk,k])
    drugGsea <-cbind(drugGsea, pathways$objects$GSEA)
    drugFet <- cbind(drugFet, pathways$objects$FET)
  }
  
  rownames(drugGsea)<-allPathways
  colnames(drugGsea)<-drugName
  
  rownames(drugFet)<-allPathways
  colnames(drugFet)<-drugName
  
  figurename1 <- paste("Hierarchical_Clustering_with_GSEA_approach_in_",modelname[k],".png",sep="")
  figurename2 <- paste("Hierarchical_Clustering_with_FET_approach_in_",modelname[k],".png",sep="")
  
  GSEA<-drugGsea
  FET<-drugFet
  
  GSEA[which(drugGsea == 0)]<-10^-4 # to avoid log(0) : 0 came from np=1000 in gsea.
  GSEA[which(drugGsea > 0.05)]<-1 # to avoid log(0) : 0 came from np=1000 in gsea.
  
  gsea<-unique(GSEA)
  
  FET[which(drugFet == 0)]<-10^-4 # to avoid log(0) : 0 came from np=1000 in gsea.
  FET[which(drugFet > 0.05)]<-1 # to avoid log(0) : 0 came from np=1000 in gsea.
  
  
  gsea<-GSEA[which(apply(GSEA,1,var) !=0),which(apply(GSEA,2,var) !=0)]
  fet<-FET[which(apply(FET,1,var) !=0),which(apply(FET,2,var) !=0)]
  
    
  #heatmap(-log10(gsea),col =  gray.colors(50, start = 1, end = 0.2, gamma = 4))
  
  
  png(figurename1,width =960, height = 960)
  GSEA_hclust<-heatmap(-log10(gsea),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10))
  dev.off()
  
  png(figurename2,width =960, height = 960)
  FET_hclust<-heatmap(-log10(fet),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10))
  dev.off()
  
  hcResult[[k]]<-list(GSEA = GSEA_hclust, FET=FET_hclust )
  allDrugStatistics[[k]]<-list(GSEA = gsea, FET=fet)
  
  aPlot <- Media(list(name = paste("GSEA_",modelname[k],sep=""), type="M", parentId=datasetID))
  aPlot <- addFile(aPlot,figurename1)
  aPlot <- storeEntity(aPlot)
  
  bPlot <- Media(list(name = paste("FET_",modelname[k],sep=""), type="M", parentId=datasetID))
  bPlot <- addFile(bPlot,figurename2)
  bPlot <- storeEntity(bPlot)
}

Result<-list(hcResult = hcResult,allDrugStatistics=allDrugStatistics)

submittedModelDatasetId <- datasetID
submittedModelLayer <- Layer(list(name = "AllStatResult", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- Layer(list(name = "AllStatResult", type = "E", parentId = submittedModelDatasetId))

submittedModelLayer <- addObject(entity=submittedModelLayer, object=Result)
submittedModelLayer <- storeEntity(submittedModelLayer)
ids<-c(ids,propertyValue(submittedModelLayer, "id"))