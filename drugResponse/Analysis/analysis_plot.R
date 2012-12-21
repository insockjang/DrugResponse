library(affy)
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF

NAME<-colnames(pData(adf_drug))

ENet<-list()
Lasso<-list()
Ridge<-list()
Svm <-list()
Pcr <-list()
Pls <-list()
Spls <-list()
for(k in 1:24){
  filenames1 <- paste("ENet_original/cvDrug_",k,"_ENet.Rdata",sep="")
  load(filenames1)
  TrPr <-c()
  TrOb <-c()
  TePr <-c()
  TeOb <-c()  
  for (k1 in 1:5){
    TrPr  <-c(TrPr,resultsOriginal[[k1]]$trainPredictions)
    TrOb  <-c(TrOb,resultsOriginal[[k1]]$trainObservations)
    TePr  <-c(TePr,resultsOriginal[[k1]]$testPredictions)
    TeOb  <-c(TeOb,resultsOriginal[[k1]]$testObservations)
  }
  ENet[[k]]<-list(trainPredictions =TrPr,trainObservations=TrOb, testPredictions = TePr, testObservations = TeOb)
  
  
  filenames2 <- paste("Lasso_original/cvDrug_",k,"_Lasso.Rdata",sep="")
  load(filenames2)
  TrPr <-c()
  TrOb <-c()
  TePr <-c()
  TeOb <-c()  
  for (k1 in 1:5){
    TrPr  <-c(TrPr,resultsOriginal[[k1]]$trainPredictions)
    TrOb  <-c(TrOb,resultsOriginal[[k1]]$trainObservations)
    TePr  <-c(TePr,resultsOriginal[[k1]]$testPredictions)
    TeOb  <-c(TeOb,resultsOriginal[[k1]]$testObservations)
  }
  Lasso[[k]]<-list(trainPredictions =TrPr,trainObservations=TrOb, testPredictions = TePr, testObservations = TeOb)
  
  filenames3 <- paste("Ridge_original/cvDrug_",k,"_Ridge.Rdata",sep="")
  load(filenames3)
  TrPr <-c()
  TrOb <-c()
  TePr <-c()
  TeOb <-c()  
  for (k1 in 1:5){
    TrPr  <-c(TrPr,resultsOriginal[[k1]]$trainPredictions)
    TrOb  <-c(TrOb,resultsOriginal[[k1]]$trainObservations)
    TePr  <-c(TePr,resultsOriginal[[k1]]$testPredictions)
    TeOb  <-c(TeOb,resultsOriginal[[k1]]$testObservations)
  }
  Ridge[[k]]<-list(trainPredictions =TrPr,trainObservations=TrOb, testPredictions = TePr, testObservations = TeOb)

  filenames4 <- paste("SVM/cvDrug_",k,"_SVM.Rdata",sep="")
  load(filenames4)
  Svm[[k]] <-list(trainPredictions =resultsScale$trainPredictions,trainObservations=resultsScale$trainObservations, 
                  testPredictions = resultsScale$testPredictions, testObservations = resultsScale$testObservations)
  
  filenames5 <- paste("PCR/cvDrug_",k,"_PCR_scale.Rdata",sep="")
  load(filenames5)
  Pcr[[k]] <-list(trainPredictions =resultsScale$trainPredictions,trainObservations=resultsScale$trainObservations, 
                testPredictions = resultsScale$testPredictions, testObservations = resultsScale$testObservations)

  filenames6 <- paste("PLS/cvDrug_",k,"_PLS_scale.Rdata",sep="")
  load(filenames6)
  Pls[[k]] <-list(trainPredictions =resultsScale$trainPredictions,trainObservations=resultsScale$trainObservations, 
                testPredictions = resultsScale$testPredictions, testObservations = resultsScale$testObservations)

  filenames7 <- paste("Spls/cvDrug_",k,"_Spls.Rdata",sep="")
  load(filenames7)
  TrPr <-c()
  TrOb <-c()
  TePr <-c()
  TeOb <-c()  
  for (k1 in 1:5){
    TrPr  <-c(TrPr,resultsOriginal[[k1]]$trainPredictions)
    TrOb  <-c(TrOb,resultsOriginal[[k1]]$trainObservations)
    TePr  <-c(TePr,resultsOriginal[[k1]]$testPredictions)
    TeOb  <-c(TeOb,resultsOriginal[[k1]]$testObservations)
  }
  Spls[[k]]<-list(trainPredictions =TrPr,trainObservations=TrOb, testPredictions = TePr, testObservations = TeOb)
}

CVresults<-list(Lasso = Lasso,ENet = ENet, Ridge = Ridge, PCR = Pcr, PLS = Pls, SVM = Svm, Spls = Spls  )



########################################################### very part
#par(mfrow = c(2,1))

C1Train<-c()
C1Test<-c()

C2Train<-c()
C2Test<-c()

C3Train<-c()
C3Test<-c()

C4Train<-c()
C4Test<-c()

C5Train<-c()
C5Test<-c()

C6Train<-c()
C6Test<-c()

C7Train<-c()
C7Test<-c()

for(k in 1:24){
  C1Train <- c(C1Train, cor(ENet[[k]]$trainObservations,ENet[[k]]$trainPredictions)^2)
  C1Test  <- c(C1Test, cor(ENet[[k]]$testObservations,ENet[[k]]$testPredictions)^2)
  
  C2Train <- c(C2Train, cor(Lasso[[k]]$trainObservations,Lasso[[k]]$trainPredictions)^2)
  C2Test  <- c(C2Test, cor(Lasso[[k]]$testObservations,Lasso[[k]]$testPredictions)^2)
  
  C3Train <- c(C3Train, cor(Ridge[[k]]$trainObservations,Ridge[[k]]$trainPredictions)^2)
  C3Test  <- c(C3Test, cor(Ridge[[k]]$testObservations,Ridge[[k]]$testPredictions)^2)
  
  C4Train <- c(C4Train, cor(Pcr[[k]]$trainObservations,Pcr[[k]]$trainPredictions)^2)
  C4Test  <- c(C4Test, cor(Pcr[[k]]$testObservations,Pcr[[k]]$testPredictions)^2)
  
  C5Train <- c(C5Train, cor(Pls[[k]]$trainObservations,Pls[[k]]$trainPredictions)^2)
  C5Test  <- c(C5Test, cor(Pls[[k]]$testObservations,Pls[[k]]$testPredictions)^2)
  
  C6Train <- c(C6Train, cor(Svm[[k]]$trainObservations,Svm[[k]]$trainPredictions)^2)
  C6Test  <- c(C6Test, cor(Svm[[k]]$testObservations,Svm[[k]]$testPredictions)^2)
  
  C7Train <- c(C7Train, cor(Spls[[k]]$trainObservations,Spls[[k]]$trainPredictions)^2)
  C7Test  <- c(C7Test, cor(Spls[[k]]$testObservations,Spls[[k]]$testPredictions)^2)

}

namesLegend<-c("Lasso","ENet","Ridge","PCR","PLS","SVM","Spls")

plot_colors <- c("black","blue","red","forestgreen","cyan","magenta","green")

plot(C1Train,col = plot_colors[1],type = "b",ylim = c(0,1),pch = 1,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = .9)

axis(1, at=1:24, lab=NAME,las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))


points(C2Train,col = plot_colors[2],type = "b",pch = 2,lty =1,cex = .9)
points(C3Train,col = plot_colors[3],type = "b",pch = 3,lty =1,cex = .9)
points(C4Train,col = plot_colors[4],type = "b",pch = 4,lty =1,cex = .9)
points(C5Train,col = plot_colors[5],type = "b",pch = 5,lty =1,cex = .9)
points(C6Train,col = plot_colors[6],type = "b",pch = 6,lty =1,cex = .9)
points(C7Train,col = plot_colors[7],type = "b",pch = 6,lty =1,cex = .9)

legend(1,.3,namesLegend,pch = seq(1,length(namesLegend)),col = plot_colors[1:length(namesLegend)],lty =1,title = "Train",cex = .7)

points(C1Test,col = plot_colors[1],type = "b",pch=1,lty =2,cex = .9)
points(C2Test,col = plot_colors[2],type = "b",pch=2,lty =2,cex = .9)
points(C3Test,col = plot_colors[3],type = "b",pch=3,lty =2,cex = .9)
points(C4Test,col = plot_colors[4],type = "b",pch=4,lty =2,cex = .9)
points(C5Test,col = plot_colors[5],type = "b",pch=5,lty =2,cex = .9)
points(C6Test,col = plot_colors[6],type = "b",pch=6,lty =2,cex = .9)
points(C7Test,col = plot_colors[7],type = "b",pch=6,lty =2,cex = .9)
legend(21,.3,namesLegend,pch = seq(1,length(namesLegend)),col = plot_colors[1:length(namesLegend)],lty =2,title = "Test",cex = .7)

title("Performance Comparison with correlation", ylab="Correlation")

analysis <- Dataset(list(
  name="Predictive Performance - How precise each data with 5 fold CV",
  description="7 predictive models",
  parentId="161936"
  ))
(dataset <- createEntity(analysis))

datasetID<-propertyValue(dataset, "id")


submittedModelDatasetId <- datasetID
submittedModelLayer <- Layer(list(name = "AllModelPerformance", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- addObject(entity=submittedModelLayer, object=CVresults)
submittedModelLayer <- storeEntity(submittedModelLayer)
