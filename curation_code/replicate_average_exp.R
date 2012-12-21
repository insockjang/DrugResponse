# Go to the directory which has raw (CEL) data you want to curate
# then, run below function
###
### library(affy)
### fit.RMA <- justRMA()
### fits$hthgu133a <- justGCMRA()


# annotation file (CEL filenames <-> SampleName_PrimarySite infos)
load("sample_annot.RData")

# find samples which have replicates
dd_annot <- unique(sample_annot[duplicated(sample_annot)])

## check which samples are replicates
annot <-unique(sample_annot)
C<-c()
for (i in 1:length(sample_annot)){
  a<-length(which(annot[i]==sample_annot))
  if (a>1) C<-rbind(C,i)
}
# make an unique sample expression by geometric mean because data already log2 transformed   
A<-c()
for (i in 1:length(dd_annot)){
  a<-which(dd_annot[i]==sample_annot)
  
  if (length(a)==2) {
    A<-cbind(A,sqrt(exprs(fits$hthgu133a)[,a[1]]*exprs(fits$hthgu133a)[,a[2]]))
    print(i)
  }
  if (length(a)==3) A<-cbind(A,(exprs(fits$hthgu133a)[,a[1]]*exprs(fits$hthgu133a)[,a[2]]*exprs(fits$hthgu133a)[,a[3]])^(1/3))
}
colnames(A)<-dd_annot

# make all unique sample (replicates are condensed by geometric mean)
d_annot <-unique(sample_annot)
single_samples_annot<-setdiff(d_annot,dd_annot)
MATRIX <-c()
for (i in 1:length(d_annot)){
  if (sum(d_annot[i] == dd_annot)) {MATRIX<-cbind(MATRIX,A[,which(d_annot[i]==dd_annot)])}
  else 
    {MATRIX <- cbind(MATRIX,exprs(fits$hthgu133a)[,which(d_annot[i]==sample_annot)])}
}
colnames(MATRIX) <- d_annot

fit <- new("ExpressionSet",exprs = MATRIX)

save(fit,file = "Sanger_exp_metaGEO.ROBJECT")

library(synapseClient)
library(affy)
synapseLogin()

## set up a project
myName <- "Sanger_Exp_metaGEO"
projName <- sprintf("%ss Curation Project %s", myName, as.character(gsub("-",".",Sys.Date())))

myProj <- createEntity(Project(list(name=projName)))

# This is parent ID(project id)
Sanger_parentID <- "5019"

## create a dataset 
myDataset <- createEntity(Dataset(list(name="Sanger Expression metaGEO", parentId=Sanger_parentID)))

## create a layer if parent id is mapped into a dataset
#myLayer <- createEntity(Layer(list(name="Sanger Expression RMA", type="E", parentId=Sanger_parentID)))

## view the dataset on the web to add a description
onWeb(myDataset)

## refresh the local copy of myDataset
myDataset <- refreshEntity(myDataset)


################
## Section 2: Working with data
################

## create a new expression layer
myExpr <- createEntity(Layer(list(name="curated expression", type="E", parentId = propertyValue(myDataset, "id"), status="curated")))

## add an annotation specifying the data format
annotValue(myExpr, "format") <- "sageBioCurated_withROBJECT"

## add the pm data file to the entity
myExpr <- addFile(myExpr, "Sanger_exp_metaGEO.ROBJECT")

## store the data
myExpr <- storeEntity(myExpr)
