# ENet<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/ENet",KK)
# Ridge<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge",KK)
# Lasso<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Lasso",KK)
# PCR<-<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PCR",KK)
# PLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PLS",KK)
# SVM<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SVM",KK)
# SPLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SPLS",KK)
source('~/COMPBIO/trunk/users/jang/drugResponse/Analysis/cvPredictivePerformance.R')
KK=24
library(affy)
library(foreach)
library(affy)
library(synapseClient)
id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))
methods= c("ENet","Ridge","Lasso","PCR","PLS","SVM","SPLS")
cvResult<-list()

cvResult[[1]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/ENet",KK)
cvResult[[2]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge",KK)
cvResult[[3]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Lasso",KK)
cvResult[[4]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PCR",KK)
cvResult[[5]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PLS",KK)
cvResult[[6]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SVM",KK)
cvResult[[7]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SPLS",KK)

a<-sort(cvResult[[1]]$overallTestCor,index.return=T)
cvMatrix<-matrix(0,nrow=5,ncol=7*KK)
a1<-seq(1,7*KK,by=7)
a2<-seq(2,7*KK,by=7)
a3<-seq(3,7*KK,by=7)
a4<-seq(4,7*KK,by=7)
a5<-seq(5,7*KK,by=7)
a6<-seq(6,7*KK,by=7)
a7<-seq(7,7*KK,by=7)

cvMatrix[,a1]<-cvResult[[1]]$indTestCor[,a$ix]
cvMatrix[,a2]<-cvResult[[2]]$indTestCor[,a$ix]
cvMatrix[,a3]<-cvResult[[3]]$indTestCor[,a$ix]
cvMatrix[,a4]<-cvResult[[4]]$indTestCor[,a$ix]
cvMatrix[,a5]<-cvResult[[5]]$indTestCor[,a$ix]
cvMatrix[,a6]<-cvResult[[6]]$indTestCor[,a$ix]
cvMatrix[,a7]<-cvResult[[7]]$indTestCor[,a$ix]

cvVector<-matrix(0,nrow=1,ncol=7*KK)
cvVector[a1]<-cvResult[[1]]$overallTestCor[a$ix]
cvVector[a2]<-cvResult[[2]]$overallTestCor[a$ix]
cvVector[a3]<-cvResult[[3]]$overallTestCor[a$ix]
cvVector[a4]<-cvResult[[4]]$overallTestCor[a$ix]
cvVector[a5]<-cvResult[[5]]$overallTestCor[a$ix]
cvVector[a6]<-cvResult[[6]]$overallTestCor[a$ix]
cvVector[a7]<-cvResult[[7]]$overallTestCor[a$ix]


plot_colors <- c("black","blue","red","forestgreen","cyan","magenta","green")
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0,0.7),pch = 15,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = 1.5)
abline(v=seq(7.5,7*(KK+1),by=7),lty =3)
axis(1, at=seq(4,7*KK,by=7), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(7*(KK-4),.3,methods,pch = 15,col = plot_colors[1:length(methods)],title = "Models",cex = 1)
