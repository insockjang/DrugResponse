# Go to the directory which has raw (CEL) data you want to curate
# then, run below function
###
### library(affy)
### fit.RMA <- justRMA()
### fit.GCRMA <- justGCMRA()

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
  
  if (length(a)==2) A<-cbind(A,sqrt(exprs(fit.GCRMA)[,a[1]]*exprs(fit.GCRMA)[,a[2]]))
  if (length(a)==3) A<-cbind(A,(exprs(fit.GCRMA)[,a[1]]*exprs(fit.GCRMA)[,a[2]]*exprs(fit.GCRMA)[,a[3]])^(1/3))
}
colnames(A)<-dd_annot

# make all unique sample (replicates are condensed by geometric mean)
d_annot <-unique(sample_annot)
single_samples_annot<-setdiff(d_annot,dd_annot)
MATRIX <-c()
for (i in 1:length(d_annot)){
  if (sum(d_annot[i] == dd_annot)) {MATRIX<-cbind(MATRIX,A[,which(d_annot[i]==dd_annot)])}
  else 
    {MATRIX <- cbind(MATRIX,exprs(fit.GCRMA)[,which(d_annot[i]==sample_annot)])}
}
colnames(MATRIX) <- d_annot

FIT.GCRMA <- new("ExpressionSet",exprs = MATRIX)

######### DONE ############
######### Save this R object into Synapse
