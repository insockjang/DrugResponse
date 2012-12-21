# ENet<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/ENet",24)
# Ridge<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge",24)
# Lasso<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Lasso",24)
# PCR<-<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PCR",24)
# PLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PLS",24)
# SVM<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SVM",24)
# SPLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SPLS",24)
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

cvResult[[1]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/ENet",24)
cvResult[[2]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge",24)
cvResult[[3]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/Lasso",24)
cvResult[[4]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PCR",24)
cvResult[[5]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/PLS",24)
cvResult[[6]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SVM",24)
cvResult[[7]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/SPLS",24)

a<-sort(cvResult[[1]]$overallTestCor,index.return=T)
cvMatrix<-matrix(0,nrow=5,ncol=7*24)
a1<-seq(1,7*24,by=7)
a2<-seq(2,7*24,by=7)
a3<-seq(3,7*24,by=7)
a4<-seq(4,7*24,by=7)
a5<-seq(5,7*24,by=7)
a6<-seq(6,7*24,by=7)
a7<-seq(7,7*24,by=7)

cvMatrix[,a1]<-cvResult[[1]]$indTestCor[,a$ix]
cvMatrix[,a2]<-cvResult[[2]]$indTestCor[,a$ix]
cvMatrix[,a3]<-cvResult[[3]]$indTestCor[,a$ix]
cvMatrix[,a4]<-cvResult[[4]]$indTestCor[,a$ix]
cvMatrix[,a5]<-cvResult[[5]]$indTestCor[,a$ix]
cvMatrix[,a6]<-cvResult[[6]]$indTestCor[,a$ix]
cvMatrix[,a7]<-cvResult[[7]]$indTestCor[,a$ix]

cvVector<-matrix(0,nrow=1,ncol=7*24)
cvVector[a1]<-cvResult[[1]]$overallTestCor[a$ix]
cvVector[a2]<-cvResult[[2]]$overallTestCor[a$ix]
cvVector[a3]<-cvResult[[3]]$overallTestCor[a$ix]
cvVector[a4]<-cvResult[[4]]$overallTestCor[a$ix]
cvVector[a5]<-cvResult[[5]]$overallTestCor[a$ix]
cvVector[a6]<-cvResult[[6]]$overallTestCor[a$ix]
cvVector[a7]<-cvResult[[7]]$overallTestCor[a$ix]

plot(t(cvVector))


plot_colors <- c("black","blue","red","forestgreen","cyan","magenta","green")
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0,1),pch = 15,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = .9)
abline(v=seq(7.5,7*25,by=7),lty =3)
axis(1, at=seq(3,7*24,by=7), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(7*21,.3,methods,pch = seq(1,length(methods)),col = plot_colors[1:length(methods)],lty =1,title = "Train",cex = .7)
