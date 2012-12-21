require(predictiveModeling)
require(affy)
require(synapseClient)
require(RColorBrewer)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

id_drugLayer <- "269024" # CCLE drugs 24 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug
drugName<-colnames(pData(adf_drug))

pathwayAnalysis<-function(filePath,pathwayDBName,modelName,drugName){
  
  # MSigDB from synapse
  
  mSigDB_annotations <- loadEntity(105363)
  mSigDB_symbolID <- loadEntity(105350)
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID  
  
  if(is.element(pathwayDBName,"KEGG")) allPathways <- mSigDB_annotations$objects$C2$KEGG
  if(is.element(pathwayDBName,"REACTOME")) allPathways <- mSigDB_annotations$objects$C2$REACTOME
  if(is.element(pathwayDBName,"BIOCARTA")) allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  
  GSEA<-c()
  for(k in 1:length(drugName)){
    filename = paste(filePath,"/",pathwayDBName,"_",modelName,"_",k,".Rdata",sep="")
    load(filename)
    gseaP<-c()
    for(k1 in 1:length(allPathways)){
      gseaP<-c(gseaP,analyticResult[[k1]]$gseaResult$p.value)
    }      
    GSEA<-cbind(GSEA,gseaP)  
  }
  
  colnames(GSEA)<-drugName
  rownames(GSEA)<-allPathways
  return(GSEA)
}


GSEA_pathway<-pathwayAnalysis("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/prior/CCLE/bootstrapRidge2","KEGG","TF",drugName)
GSEA<-GSEA_pathway
GSEA<-unique(GSEA)
GSEA[which(GSEA == 0)]<-10^-3 # to avoid log(0) : 0 came from np=1000 in gsea.
GSEA[which(GSEA == "NaN")] <- 1 # to avoid log(0) : 0 came from np=1000 in gsea.
GSEA_hclust<-heatmap(-log10(GSEA),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10),main = "KEGG")

plot(GSEA_hclust$Colv)


GSEA_hclust<-heatmap(-log10(A),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10),main = "KEGG")


GSEA1<-GSEA_KEGG2
GSEA1[which(GSEA1 == 0)]<-10^-3 # to avoid log(0) : 0 came from np=1000 in gsea.
GSEA1[which(GSEA1 == "NaN")] <- 1 # to avoid log(0) : 0 came from np=1000 in gsea.

GSEA1<-GSEA_REACTOME
GSEA1[which(GSEA1 == 0)]<-10^-3 # to avoid log(0) : 0 came from np=1000 in gsea.
GSEA1[which(GSEA1 == "NaN")] <- 1 # to avoid log(0) : 0 came from np=1000 in gsea.

GSEA2<-GSEA_BIOCARTA
GSEA2[which(GSEA2 == 0)]<-10^-3 # to avoid log(0) : 0 came from np=1000 in gsea.
GSEA2[which(GSEA2 == "NaN")] <- 1 # to avoid log(0) : 0 came from np=1000 in gsea.

png("Sanger_KEGG.png",width =1800, height = 1800)
GSEA_hclust<-heatmap(-log10(GSEA),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10),main = "KEGG")
dev.off()

png("Sanger_REACTOME.png",width =1800, height = 1800)
GSEA1_hclust<-heatmap(-log10(GSEA1),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10),main = "REACTOME")
dev.off()

png("Sanger_BIOCARTA.png",width =1800, height = 1800)
GSEA2_hclust<-heatmap(-log10(GSEA2),keep.dendro=TRUE,col = colorRampPalette(brewer.pal(5,"Blues"))(10),main = "BIOCARTA")
dev.off()


par(mfrow=c(2,1))
plot(GSEA_hclust$Colv,main = "single_filtering")
plot(GSEA1_hclust$Colv,main = "single_nofiltering")

#heatmap(-log10(gsea),col =  gray.colors(50, start = 1, end = 0.2, gamma = 4))


png(figurename1,width =960, height = 960)
dev.off()

##################################################################################################################
# synapseID  
##################################################################################################################

analysis <- Dataset(list(
  name="New Pathway Analysis Comparison",
  description="GSEA and FET with CCLE Drug Response depending upon predictiveModel",
  parentId="161936"
  ))
(dataset <- createEntity(analysis))

datasetID<-propertyValue(dataset, "id")


aPlot <- Media(list(name = paste("GSEA_",modelname[k],sep=""), type="M", parentId=datasetID))
aPlot <- addFile(aPlot,figurename1)
aPlot <- storeEntity(aPlot)


submittedModelDatasetId <- datasetID
submittedModelLayer <- Layer(list(name = "AllStatResult", type = "E", parentId = submittedModelDatasetId))

submittedModelLayer <- addObject(entity=submittedModelLayer, object=Result)
submittedModelLayer <- storeEntity(submittedModelLayer)
