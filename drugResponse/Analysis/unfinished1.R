library(affy)
library(synapseClient)
id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

drugName<-colnames(pData(adf_drug))

setwd("~/COMPBIO/trunk/users/jang/drugResponse/CCLE/")

load("ENet/cv5performance.Rdata")
ENet<-corPearsonScale[,2]

load("Lasso/cv5performance.Rdata")
Lasso<-corPearsonScale[,2]

load("Ridge/cv5performance.Rdata")
Ridge<-corPearsonScale[,2]

load("PLS/cv5performance.Rdata")
PLS<-corPearsonScale[,2]

load("SVM/cv5performance.Rdata")
SVM<-corPearsonScale[,2]

load("SPLS/cv5performance.Rdata")
SPLS<-corPearsonScale[,2]

COR<-rbind(ENet,Lasso,Ridge,PLS,SVM,SPLS)
colnames(COR)<-drugName

a<-sort(COR[1,],index.return = T)
COR_sort<-COR[,a$ix]


namesLegend<-rownames(COR_sort)
plot_colors <- c("black","blue","red","forestgreen","cyan","magenta","green")

plot(COR_sort[1,],col = plot_colors[1],type = "b",ylim = c(0,1),pch = 1,lty =1, ylab = "Correlation", axe = FALSE, ann=FALSE,cex = .9)
axis(1, at=1:24, lab=colnames(COR_sort),las = 2, cex.axis =0.7)
axis(2,las =1, at = seq(0,1,0.2))

points(COR_sort[2,],col = plot_colors[2],type = "b",pch = 2,lty =1,cex = .9)
points(COR_sort[3,],col = plot_colors[3],type = "b",pch = 3,lty =1,cex = .9)
points(COR_sort[4,],col = plot_colors[4],type = "b",pch = 4,lty =1,cex = .9)
points(COR_sort[5,],col = plot_colors[5],type = "b",pch = 5,lty =1,cex = .9)
points(COR_sort[6,],col = plot_colors[6],type = "b",pch = 6,lty =1,cex = .9)
points(COR_sort[7,],col = plot_colors[7],type = "b",pch = 6,lty =1,cex = .9)

legend(21,.3,namesLegend,pch = seq(1,length(namesLegend)),col = plot_colors[1:length(namesLegend)],lty =1,title = "Train",cex = .7)


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



