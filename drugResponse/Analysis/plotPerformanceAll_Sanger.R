# ENet<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/ENet",KK)
# Ridge<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Ridge",KK)
# Lasso<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Lasso",KK)
# PCR<-<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PCR",KK)
# PLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PLS",KK)
# SVM<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SVM",KK)
# SPLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SPLS",KK)
source('~/COMPBIO/trunk/users/jang/drugResponse/Analysis/cvPredictivePerformance.R')
KK=130
K=6
library(affy)
library(foreach)
library(affy)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
id_drugLayer <- "220680" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))
methods= c("ENet","Ridge","Lasso","PCR","PLS","SVM")
cvResult<-list()

cvResult[[1]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/ENet2",KK)
cvResult[[2]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Ridge2",KK)
cvResult[[3]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Lasso2",KK)
cvResult[[4]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PCR",KK)
cvResult[[5]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PLS",KK)
cvResult[[6]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SVM",KK)
#cvResult[[7]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SPLS",KK)

a<-sort(cvResult[[1]]$overallTestCor,index.return=T)
cvMatrix<-matrix(0,nrow=5,ncol=K*KK)
a1<-seq(1,K*KK,by=K)
a2<-seq(2,K*KK,by=K)
a3<-seq(3,K*KK,by=K)
a4<-seq(4,K*KK,by=K)
a5<-seq(5,K*KK,by=K)
a6<-seq(6,K*KK,by=K)
#a7<-seq(7,7*KK,by=7)

cvMatrix[,a1]<-cvResult[[1]]$indTestCor[,a$ix]
cvMatrix[,a2]<-cvResult[[2]]$indTestCor[,a$ix]
cvMatrix[,a3]<-cvResult[[3]]$indTestCor[,a$ix]
cvMatrix[,a4]<-cvResult[[4]]$indTestCor[,a$ix]
cvMatrix[,a5]<-cvResult[[5]]$indTestCor[,a$ix]
cvMatrix[,a6]<-cvResult[[6]]$indTestCor[,a$ix]
#cvMatrix[,a7]<-cvResult[[7]]$indTestCor[,a$ix]

cvVector<-matrix(0,nrow=1,ncol=K*KK)
cvVector[a1]<-cvResult[[1]]$overallTestCor[a$ix]
cvVector[a2]<-cvResult[[2]]$overallTestCor[a$ix]
cvVector[a3]<-cvResult[[3]]$overallTestCor[a$ix]
cvVector[a4]<-cvResult[[4]]$overallTestCor[a$ix]
cvVector[a5]<-cvResult[[5]]$overallTestCor[a$ix]
cvVector[a6]<-cvResult[[6]]$overallTestCor[a$ix]
#cvVector[a7]<-cvResult[[7]]$overallTestCor[a$ix]



plot_colors <- c("black","blue","red","forestgreen","cyan","magenta")
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0,0.7),pch = 15,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = 2)
abline(v=seq(7.5,K*(KK+1),by=K),lty =3)
axis(1, at=seq(3,K*KK,by=K), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
#legend(K*(KK-15),.2,methods[1:K],pch = 15,text.col = plot_colors[1:K],title = "Models",cex = 1.5)
#legend(locatK*(KK-20),.35,methods[1:K],text.col = plot_colors[1:K],cex = 1.5,bg="transparent",bty="n")
legend("bottomright",methods[1:K],text.col = plot_colors[1:K],cex = 2.0,bg="transparent",bty="n")
#legend(7*(KK-15),.3,methods,pch = seq(1,length(methods)),col = plot_colors[1:length(methods)],lty =1,title = "Train",cex = .7)
