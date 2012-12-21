setwd("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/Sanger/")

par(mfrow=c(1,2))
load("Lasso/cv5performance.Rdata")
a1<-corPearsonScale[,2]
plot(corPearsonScale[,2])
load("NP/NPLasso/cv5performance.Rdata")
a2<-corPearsonScale[,2]
points(corPearsonScale[,2],col="red")

load("Ridge/cv5performance.Rdata")
b1<-corPearsonScale[,2]
plot(corPearsonScale[,2])
load("NP/NPRidge/cv5performance.Rdata")
b2<-corPearsonScale[,2]
points(corPearsonScale[,2],col="red")

load("ENet/cv5performance.Rdata")
c1<-corPearsonScale[,2]
plot(corPearsonScale[,2])
load("NP/NPENet/cv5performance.Rdata")
c2<-corPearsonScale[,2]
points(corPearsonScale[,2],col="red")


load("SVM/cv5performance.Rdata")
d1<-corPearsonScale[,2]
plot(corPearsonScale[,2])
load("NP/NPRidge/cv5performance.Rdata")
d2<-corPearsonScale[,2]
points(corPearsonScale[,2],col="red")
dev.off()

plot(a1)
points(b1,col="red")
points(c1,col="blue")
points(d1,col="green")

xx<-sort(b1,index.return=T)
plot(b1[xx$ix])
points(a1[xx$ix],col = "red")
points(c1[xx$ix],col = "blue")
points(d1[xx$ix],col = "green")
legend(1,1,"Ridge","Lasso","ENet","SVM")

################################################################## bs Lasso
sumBSLasso<-list()
sumBSRidge<-list()
sumBSENet<-list()

for(kk in 1:KK){
  fileLasso = paste("bootstrapLasso/cvDrug_",kk,".Rdata",sep="")
  
  load(fileLasso)
  
  # Lasso processing to find bootstrapping feature selection
  ResultLasso<-c()
  for(k in 1:length(resultsLasso)){
    ResultLasso<-cbind(ResultLasso,as.matrix(resultsLasso[[k]]))
  }
  
  binResultLasso <- ResultLasso!=0
  sumBSLasso[[k]] <- apply(binResultLasso,1,sum)
}

################################################################## bs Ridge
for(kk in 1:24){
  fileRidge = paste("bootstrapRidge/cvDrug_",kk,".Rdata",sep="")
  
  load(fileRidge)
  
  # Lasso processing to find bootstrapping feature selection
  ResultRidge<-c()
  for(k in 1:length(resultsRidge)){
    B<-rank(abs(as.matrix(resultsRidge[[k]])))
    b<-B/9326
    ResultRidge<-cbind(ResultRidge,b)
  }
  
  sumBSRidge[[kk]] <- apply(ResultRidge,1,sum)
}

################################################################## bs ENet
for(kk in 1:24){
  fileENet = paste("bootstrapENet/cvDrug_",kk,".Rdata",sep="")
  
  load(fileENet)
  
  # Lasso processing to find bootstrapping feature selection
  ResultENet<-c()
  for(k in 1:length(resultsENet)){
    ResultENet<-cbind(ResultENet,as.matrix(resultsENet[[k]]))
  }
  
  binResultENet <- ResultENet!=0
  sumBSENet[[kk]] <- apply(binResultENet,1,sum)
}

