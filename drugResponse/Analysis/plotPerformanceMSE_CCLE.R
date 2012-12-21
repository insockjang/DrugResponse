library(predictiveModeling)
library(synapseClient)
library(affy)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

# CCLE
id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))

mse<-function(resultsScale){
  nfolds <- length(resultsScale)
  trPred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainPredictions}
  tePred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testPredictions}
  trObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainObservations}
  teObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  return(mean((allTeObsr - allTePred)^2))
}

correlation<-function(resultsScale){
  nfolds <- length(resultsScale)
  trPred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainPredictions}
  tePred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testPredictions}
  trObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainObservations}
  teObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  return(cor(allTeObsr,allTePred))
}
modelNum =7

MSE<-matrix(0,nrow = length(drugName), ncol = modelNum)
COR<-matrix(0,nrow = length(drugName), ncol = modelNum)
for(kk in 1:24){
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/PCR/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,1]<-mse(resultsScale)
  COR[kk,1]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/PLS/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,2]<-mse(resultsScale)
  COR[kk,2]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/SVM/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,3]<-mse(resultsScale)
  COR[kk,3]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/SPLS/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,4]<-mse(resultsScale)
  COR[kk,4]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/Lasso2/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,5]<-mse(resultsScale)
  COR[kk,5]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge2/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,6]<-mse(resultsScale)
  COR[kk,6]<-correlation(resultsScale)
  
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/ENet2/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  MSE[kk,7]<-mse(resultsScale)
  COR[kk,7]<-correlation(resultsScale)
  
}

colnames(MSE)<-c("PCR","PLS","SVM","SPLS","Lasso","Ridge","ENet")
rownames(MSE)<-drugName
colnames(COR)<-c("PCR","PLS","SVM","SPLS","Lasso","Ridge","ENet")
rownames(COR)<-drugName


KK<-nrow(MSE)
a<-sort(MSE[,7],decreasing =T,index.return=T)

a1<-seq(1,modelNum*KK,by=modelNum)
a2<-seq(2,modelNum*KK,by=modelNum)
a3<-seq(3,modelNum*KK,by=modelNum)
a4<-seq(4,modelNum*KK,by=modelNum)
a5<-seq(5,modelNum*KK,by=modelNum)
a6<-seq(6,modelNum*KK,by=modelNum)
a7<-seq(7,modelNum*KK,by=modelNum)

cvVector<-matrix(0,nrow=1,ncol=modelNum*KK)
cvVector[a1]<-MSE[a$ix,1]
cvVector[a2]<-MSE[a$ix,2]
cvVector[a3]<-MSE[a$ix,3]
cvVector[a4]<-MSE[a$ix,4]
cvVector[a5]<-MSE[a$ix,5]
cvVector[a6]<-MSE[a$ix,6]
cvVector[a7]<-MSE[a$ix,7]
png("predictivePerformanceAnalysis/CCLE_MSE.png",width = 900, height = 900)
plot_colors <- c("black","blue","red","forestgreen","magenta","cyan","green")
methods = colnames(MSE)
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0.5,1),pch = 15,lty =1, ylab = "MSE", axe = FALSE, ann=FALSE,cex = 1.5)
abline(v=seq(modelNum+0.5,modelNum*(KK+1),by=modelNum),lty =3)
axis(1, at=seq(round(modelNum/2)+0.5,modelNum*KK,by=modelNum), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(modelNum,0.73,methods,pch = 15,col = plot_colors,title = "Models",cex = 1)
title("CCLE MSE")
dev.off()

a<-sort(COR[,7],index.return=T)

a1<-seq(1,modelNum*KK,by=modelNum)
a2<-seq(2,modelNum*KK,by=modelNum)
a3<-seq(3,modelNum*KK,by=modelNum)
a4<-seq(4,modelNum*KK,by=modelNum)
a5<-seq(5,modelNum*KK,by=modelNum)
a6<-seq(6,modelNum*KK,by=modelNum)
a7<-seq(7,modelNum*KK,by=modelNum)

cvVector<-matrix(0,nrow=1,ncol=modelNum*KK)
cvVector[a1]<-COR[a$ix,1]
cvVector[a2]<-COR[a$ix,2]
cvVector[a3]<-COR[a$ix,3]
cvVector[a4]<-COR[a$ix,4]
cvVector[a5]<-COR[a$ix,5]
cvVector[a6]<-COR[a$ix,6]
cvVector[a7]<-COR[a$ix,7]
png("predictivePerformanceAnalysis/CCLE_COR.png",width = 900, height = 900)
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0.2,0.7),pch = 15,lty =1, ylab = "MSE", axe = FALSE, ann=FALSE,cex = 1.5)
abline(v=seq(modelNum+0.5,modelNum*(KK+1),by=modelNum),lty =3)
axis(1, at=seq(round(modelNum/2)+0.5,modelNum*KK,by=modelNum), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(modelNum,0.65,methods,pch = 15,col = plot_colors,title = "Models",cex = 1)
title("CCLE Correlation")
dev.off()