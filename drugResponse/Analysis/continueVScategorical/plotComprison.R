library(predictiveModeling)
library(synapseClient)
library(affy)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))

AUC<-matrix(0,nrow = 24,ncol=3)
for(kk in 1:24){
  filename1 <-paste("~/PredictiveModel/drugResponse/CCLE/catLasso2/ROCR_mean_",kk,".Rdata",sep ="")
  load(filename1)
  AUC[kk,1]<-resultsCat$erAUC@y.values[[1]]
  filename2 <-paste("~/PredictiveModel/drugResponse_woExp/CCLE/catLasso2/ROCR_mean_",kk,".Rdata",sep ="")
  load(filename2)
  AUC[kk,2]<-resultsCat$erAUC@y.values[[1]]
  filename3 <-paste("~/PredictiveModel/drugResponse_woCopy/CCLE/catLasso2/ROCR_mean_",kk,".Rdata",sep ="")
  load(filename3)  
  AUC[kk,3]<-resultsCat$erAUC@y.values[[1]]  
}

colnames(AUC)<-c("EXP+CNV+MUT","CNV+MUT","EXP+MUT")
rownames(AUC)<-drugName


a<-sort(AUC[,1],index.return=T)

# par(mfrow=c(1,3))
plot(AUC[a$ix,1],col="blue",pch = 15,lty =1,ylim = c(0.4,1),axe = FALSE, ann=FALSE,cex = 1.5)
points(AUC[a$ix,2],col="red",cex = 1.5,pch = 15,lty =1)
points(AUC[a$ix,3],col="green",cex = 1.5,pch = 15,lty =1)
axis(1, at=seq(1,length(drugName)), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(.5,1,colnames(AUC),pch = 15,col = c("blue","red","green"),title = "Models",cex = 1)
title("Categorical Lasso")
