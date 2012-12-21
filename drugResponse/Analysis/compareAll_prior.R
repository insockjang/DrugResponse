# ENet<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/ENet",KK)
# Ridge<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Ridge",KK)
# Lasso<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Lasso",KK)
# PCR<-<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PCR",KK)
# PLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/PLS",KK)
# SVM<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SVM",KK)
# SPLS<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/SPLS",KK)
KK=130
library(affy)
library(foreach)
library(affy)
library(synapseClient)
id_drugLayer <- "220680" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))
methods= c("ENet","Ridge","Lasso","ENetPrior","RidgePrior","LassoPrior")
cvResult<-list()

cvResult[[1]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/ENet",KK)
cvResult[[4]]<-cvPredictivePerformancePrior("~/COMPBIO/trunk/users/jang/drugResponse/prior/Sanger/ENet",KK)

cvResult[[2]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Ridge",KK)
cvResult[[5]]<-cvPredictivePerformancePrior("~/COMPBIO/trunk/users/jang/drugResponse/prior/Sanger/Ridge",KK)

cvResult[[3]]<-cvPredictivePerformance("~/COMPBIO/trunk/users/jang/drugResponse/Sanger/Lasso",KK)
cvResult[[6]]<-cvPredictivePerformancePrior("~/COMPBIO/trunk/users/jang/drugResponse/prior/Sanger/Lasso",KK)

a<-sort(cvResult[[1]]$overallTestCor,index.return=T)
b<-sort(cvResult[[5]]$overallTestCor,index.return=T)

par(mfrow=c(1,2))
plot(cvResult[[1]]$overallTestCor[a$ix],col="blue")
points(cvResult[[2]]$overallTestCor[a$ix],col="red")
points(cvResult[[3]]$overallTestCor[a$ix],col="green")

plot(cvResult[[4]]$overallTestCor[b$ix],col="blue")
plot(cvResult[[5]]$overallTestCor[b$ix],col="red")
points(cvResult[[6]]$overallTestCor[b$ix],col="green")



plot_colors <- c("black","blue","red","forestgreen","cyan","magenta")
plot(t(cvVector),col = plot_colors,type = "p",ylim = c(0,0.7),pch = 15,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = .9)
abline(v=seq(6.5,6*(KK+1),by=6),lty =3)
axis(1, at=seq(3,6*KK,by=6), lab=drugName[a$ix],las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))
legend(7*(KK-15),.3,methods,pch = 15,col = plot_colors[1:length(methods)],title = "Models",cex = 1)
#legend(7*(KK-15),.3,methods,pch = seq(1,length(methods)),col = plot_colors[1:length(methods)],lty =1,title = "Train",cex = .7)
