TF<-read.delim("tf_list_nov08_entrez_1875.txt",header = F)
M<-read.table("affy_sglist_u133p2_na27.txt",header = F)

library(hgu133plus2.db)

# Symbol and probeset mapping ###########################
x <- hgu133plus2SYMBOL
y <- hgu133plus2ENTREZID
# Get the probe identifiers that are mapped to a gene symbol

xx<-list()
yy<-list()
for(k in 1:nrow(M)){
# mapped_probe <- mappedkeys(as.character(M[k,1]))
# Convert to a list
  xx[[k]] <- as.list(x[as.character(M[k,1])])
  yy[[k]] <- as.list(y[as.character(M[k,1])])
}

probeGeneMap <-c()
for(i in 1:length(xx)){
  probeGeneMap <- rbind(probeGeneMap,c(as.character(names(xx[[i]])),as.character(xx[[i]]),as.character(yy[[i]])))
}
write.table(probeGeneMap[,c(3,2,1)],"modulator_list.txt",row.names = F, col.names=F)
