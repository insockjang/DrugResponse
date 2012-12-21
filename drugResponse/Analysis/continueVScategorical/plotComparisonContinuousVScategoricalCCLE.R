library(predictiveModeling)
library(synapseClient)
library(affy)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))

AUC<-matrix(0,nrow = 24,ncol=8)
for(kk in 1:24){
  filename1 <-paste("~/PredictiveModel/drugResponse/CCLE/PCR/ROCR_",kk,".Rdata",sep ="")
  load(filename1)
  AUC[kk,1]<-resultsCat$erAUC@y.values[[1]]
  filename2 <-paste("~/PredictiveModel/drugResponse/CCLE/PLS/ROCR_",kk,".Rdata",sep ="")
  load(filename2)
  AUC[kk,2]<-resultsCat$erAUC@y.values[[1]]
  filename3 <-paste("~/PredictiveModel/drugResponse/CCLE/SVM/ROCR_",kk,".Rdata",sep ="")
  load(filename3)
  AUC[kk,3]<-resultsCat$erAUC@y.values[[1]]
  filename4 <-paste("~/PredictiveModel/drugResponse/CCLE/SPLS/ROCR_",kk,".Rdata",sep ="")
  load(filename4)
  AUC[kk,4]<-resultsCat$erAUC@y.values[[1]]
  
  filename5 <-paste("~/PredictiveModel/drugResponse/CCLE/Lasso2/ROCR_",kk,".Rdata",sep ="")
  load(filename5)  
  AUC[kk,5]<-resultsCat$erAUC@y.values[[1]]  
  filename6 <-paste("~/PredictiveModel/drugResponse/CCLE/Ridge2/ROCR_",kk,".Rdata",sep ="")
  load(filename6)
  AUC[kk,6]<-resultsCat$erAUC@y.values[[1]]
  filename7 <-paste("~/PredictiveModel/drugResponse/CCLE/ENet2/ROCR_",kk,".Rdata",sep ="")
  load(filename7)
  AUC[kk,7]<-resultsCat$erAUC@y.values[[1]]
  filename8 <-paste("~/PredictiveModel/drugResponse/CCLE/randomForest/ntree100/ROCR_",kk,".Rdata",sep ="")
  load(filename8)
  AUC[kk,8]<-resultsCat$erAUC@y.values[[1]]
  
}

colnames(AUC)<-c("PCR","PLS","SVM","SPLS","Lasso","Ridge","ENet","categoricalRandomForest")
rownames(AUC)<-drugName


modelNum = 8
KK<-nrow(AUC)
a<-sort(AUC[,7],index.return=T)

a1<-seq(1,modelNum*KK,by=modelNum)
a2<-seq(2,modelNum*KK,by=modelNum)
a3<-seq(3,modelNum*KK,by=modelNum)
a4<-seq(4,modelNum*KK,by=modelNum)
a5<-seq(5,modelNum*KK,by=modelNum)
a6<-seq(6,modelNum*KK,by=modelNum)
a7<-seq(7,modelNum*KK,by=modelNum)
a8<-seq(8,modelNum*KK,by=modelNum)

cvVector<-matrix(0,nrow=1,ncol=modelNum*KK)
cvVector[a1]<-AUC[a$ix,1]
cvVector[a2]<-AUC[a$ix,2]
cvVector[a3]<-AUC[a$ix,3]
cvVector[a4]<-AUC[a$ix,4]
cvVector[a5]<-AUC[a$ix,5]
cvVector[a6]<-AUC[a$ix,6]
cvVector[a7]<-AUC[a$ix,7]
cvVector[a8]<-AUC[a$ix,8]

png("predictivePerformanceAnalysis/CCLE_AUC.png",width = 900, height = 900)
plot_colors <- c("black","blue","red","forestgreen","magenta","cyan","green","darkred")
methods = colnames(AUC)
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0.5,0.85),pch = 15,lty =1, ylab = "AUC", axe = FALSE, ann=FALSE,cex = 1.5)
abline(v=seq(modelNum+0.5,modelNum*(KK+1),by=modelNum),lty =3)
axis(1, at=seq(round(modelNum/2)+0.5,modelNum*KK,by=modelNum), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(modelNum,0.85,methods,pch = 15,col = plot_colors,title = "Models",cex = 1)
title("CCLE ROC AUC")
dev.off()