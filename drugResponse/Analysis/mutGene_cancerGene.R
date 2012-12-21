targets<-Compound$Targets
cancerGenes<-c()
for(k in 1:length(targets)){
  mydata1 <- gsub("[ ]","",as.character(targets[[k]]))
  mydata2 <- strsplit(mydata1,",",fixed=TRUE)[[1]]
  cancerGenes<-union(cancerGenes,mydata2)
}


id_oncomapLayer <- "266141"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap


mutGenes<-rownames(exprs(eSet_oncomap))

K <-c()
for(k in 1:47){
  aa<-grep(mutGenes[k],toupper(cancerGenes))
  if(length(aa)>0) K<-c(K,k)
}