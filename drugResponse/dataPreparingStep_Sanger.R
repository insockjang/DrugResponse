library(gdata)
file_expr<-"/home/ijang/PredictiveModel/Sanger/Sanger800_expr_rma_entrez.txt"
file_annot <- "/home/ijang/PredictiveModel/Sanger/Sanger_affy_n798_sample_info_published.xls"

annot<-read.xls(file_annot)
expr<-read.table(file_expr)

################################################ map cel file -> cell line names within annotation file
# column name is cel file name -> cell line names
sampleAnnotExpr<-sub(".CEL","",names(expr))

annotSampleName<-annot$SampleName
annotSampleName<-gsub("[-]","",annotSampleName)
annotSampleCel<-make.names(annot$Scan)

A<-match(sampleAnnotExpr,annotSampleCel) # which sampleAnnotExpr belongs to annotSampleCel
names(expr)<-annotSampleName[A]

################################################## map EntrezID into geneSymbol
geneAnnot<-rownames(expr)
geneAnnot<-sub("_at","",geneAnnot)

# org.Hs.egGENENAME Map between Entrez Gene IDs and Genes
library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x)
# Convert to a list

y<-is.element(geneAnnot,mapped_genes)
y1<-which(y!=1)

geneSymbolAnnot<-c()
for(k in 1:length(geneAnnot)){  
  if(is.element(k,y1)){
    geneSymbolAnnot<-c(geneSymbolAnnot,paste("LOC",geneAnnot[k],sep=""))
  }
  else{
    xx <- as.list(x[geneAnnot[k]])
    geneSymbolAnnot<-c(geneSymbolAnnot,as.character(xx))
  }  
}

rownames(expr)<-geneSymbolAnnot
expr2<-expr[-y1,] # Non annotated gene EntrezID should be filtered out
expr2<-as.matrix(expr2) # semifinal result but I need to check sample-replicates


################################################ replicated samples should be handled on average 
Names<-colnames(expr2)
uniqueNames<-unique(Names)
mean.gm<-function(x) exp(mean(log(x)))

uniqueExpr<-c()
for(k in 1:length(uniqueNames)){
  a<-match(uniqueNames[k],Names)
  if(length(a)==1){
    uniqueExpr<-cbind(uniqueExpr,expr2[,a])
  }
  else{
    uniqueExpr<-cbind(uniqueExpr,apply(expr2[,a],1,mean.gm)) # geometric mean is applied becasue of being already log2 transformed
  }  
}


colnames(uniqueExpr)<-uniqueNames
eSet_expr <-ExpressionSet(uniqueExpr)





########################################################################################## copy number

file_copy<-"/home/ijang/PredictiveModel/Sanger/Sanger800_copy_number_by_gene.txt"

annot<-read.xls(file_annot)
copy2<-read.delim2(file_copy)

Name<-names(copy2)[-1]
B<-strsplit(Name,"_")

BB<-c()
for(k in 1:length(B)){
  BB<-c(BB,B[[k]][1])
}
BB<-gsub("[.]","",BB)
BB<-make.names(BB)

copy<-copy2[,-1]
names(copy)<-BB
rownames(copy)<-copy2[,1]

copy1<-c()

for(k2 in 1:ncol(copy)){
  copy1<-cbind(copy1,as.numeric(as.character(copy[,k2])))  
}

rownames(copy1)<-rownames(copy)
colnames(copy1)<-names(copy)

eSet_copy <-ExpressionSet(as.matrix(copy1))
############################################################################################


############################################################################### oncomap

file_mut1 <- "/home/ijang/PredictiveModel/Sanger/TCN_cell_lines_120310.csv"

mut<-read.csv(file_mut1)

geneName<-names(mut)[-1]
cellLineName<-mut$X.gene[-c(1:3)]
cellLineName<-gsub("[-]","",cellLineName)
cellLineName<-make.names(cellLineName)

mut2<-as.matrix(mut[-c(1,2,3),-1])

rownames(mut2)<-cellLineName
colnames(mut2)<-names(mut)[-1]


mut1<-c()
for(k2 in 1:ncol(mut2)){
  mut1<-c(mut1,as.numeric(as.character(mut2[,k2])))
}

rownames(mut1)<-rownames(mut2)
colnames(mut1)<-colnames(mut2)

eSet_oncomap <-ExpressionSet(mut1)
############################################################################################## 
file_mut1 <- "/home/ijang/PredictiveModel/Sanger/CosmicCellLineExport_v58_150312.tsv"

mut<-read.delim(file_mut1)

cellLineName<-mut$Sample.name[2:2037]
cellLineName<-gsub("[-]","",cellLineName)
cellLineName<-make.names(cellLineName)

geneName<-as.character(mut$Gene.name[2:2037])

mat<-matrix(0,nrow = length(unique(geneName)),ncol=length(unique(cellLineName)))
rownames(mat)<-as.character(unique(geneName))
colnames(mat)<-as.character(unique(cellLineName))
for(k in 1:length(geneName)){
  mat[geneName[k],cellLineName[k]]<-1
}

eSet_oncomap <-ExpressionSet(mat)

############################################################################### mut

file_mut2 <- "/home/ijang/PredictiveModel/Sanger/CosmicCompleteExport_v58_150312.tsv"

mut1<-read.delim(file_mut2)

cellLineName1<-mut1$Sample.name
cellLineName1<-gsub("[-]","",cellLineName1)
cellLineName1<-make.names(cellLineName1)
geneName1<-mut1$Gene.name

cellLineName<-make.names(annotSampleName)
mat<-matrix(0,nrow = length(unique(geneName1)),ncol=length(cellLineName))
rownames(mat)<-as.character(unique(geneName1))
colnames(mat)<-as.character(cellLineName)
for(k in 1:length(geneName1)){
  if(!is.na(match(cellLineName1[k],cellLineName)))  mat[geneName1[k],cellLineName1[k]]<-1
}
a1<-apply(mat,1,sum)
A1<-which(a1==0)
a2<-apply(mat,2,sum)
A2<-which(a2==0)
mat<-mat[-A1,-A2]
eSet_mut <-ExpressionSet(mat)

##################################################################################### cancerGeneMut

file_cancer <- "/home/ijang/PredictiveModel/Sanger/cancer_gene_census.xls"
cancer<-read.xls(file_cancer)
cancerGeneName<-as.character(cancer$Symbol)

mat1<-mat[intersect(rownames(mat),cancerGeneName),]
a1<-apply(mat1,1,sum)
A1<-which(a1==0)
a2<-apply(mat1,2,sum)
A2<-which(a2==0)
mat2<-mat1[-A1,-A2]
eSet_cancerGeneMut <- ExpressionSet(mat2)


####################################################################### oncomap
file1 = "/home/ijang/PredictiveModel/Sanger/gdsc_manova_input_w1.xls"
sanger<-read.xls(file1)

cellLineName<-as.character(sanger$Cell.Line[-1])
cellLineName<-make.names(gsub("[-]","",cellLineName))

drugName<-names(sanger)

geneName<-drugName[6:73]

mut<-c()
for(k in 1:length(geneName)){
  b<-grep("p.",as.character(sanger[[5+k]]))
  d<-grep("na",as.character(sanger[[5+k]]))
  a<-rep(0,length(sanger[[6]]))
  a[b]<-1  
  a[d]<-NA
  mut<-cbind(mut,a)
}
MUT<-mut[-1,]
colnames(MUT)<-geneName
rownames(MUT)<-cellLineName
tMUT <-t(MUT)

KK<-ncol(tMUT)
K<-c()
for(k in 1:nrow(tMUT)){
  if(sum(is.na(tMUT[k,]))==KK){
    K<-c(K,k) 
  }
  else{
    if(sum(tMUT[k,])==0) K<-c(K,k)
  }
  
}

geneMUT<-tMUT[-K,]

eSet_oncomap<-ExpressionSet(geneMUT)

############################################################################### drug

file1 = "/home/ijang/PredictiveModel/Sanger/gdsc_manova_input_w1.xls"
sanger<-read.xls(file1)

cellLineName<-as.character(sanger$Cell.Line[-1])
cellLineName<-make.names(gsub("[-]","",cellLineName))

drugName<-names(sanger)


a<-grep("IC_50",drugName)

dataSanger<-as.matrix(sanger[-1,a])
drugName_ic50<-drugName[a]

B<-strsplit(drugName_ic50,"_")

drugName1<-c()
for(k in 1:length(B)){
  drugName1<-c(drugName1,B[[k]][1])
}

drugName1<-gsub("[.]","-",drugName1)[1:131]

dataSanger<-dataSanger[,1:131]

rownames(dataSanger) <- cellLineName
colnames(dataSanger) <- drugName1

dataSanger<-dataSanger[,-82] # Camptothecin is replicated but the measured values are different



dataSanger1<-c()
for(k2 in 1:ncol(dataSanger)){
  dataSanger1<-cbind(dataSanger1,as.numeric(as.character(dataSanger[,k2])))
}

rownames(dataSanger1)<-rownames(dataSanger)
colnames(dataSanger1)<-colnames(dataSanger)



adf_drug<-AnnotatedDataFrame(data=as.data.frame(dataSanger1))




dataset <- Dataset(list(
  name="Sanger",
  parentId=170143
  ))
dataset <- createEntity(dataset)

#onWeb(dataset)

probeIDLayer <- Layer(list(name = "R_Expr", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_expr)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Cnv", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_copy)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Oncomap", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_oncomap)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_CancerMut", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, eSet_cancerGeneMut)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Drug", type = "E", parentId = propertyValue(dataset, "id"), status="db"))
testLayer <- addObject(probeIDLayer, adf_drug)
testLayer <- storeEntity(testLayer)

