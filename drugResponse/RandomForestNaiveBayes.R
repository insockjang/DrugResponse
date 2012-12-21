# DEMO naive Bayes Classification

library(class)
library(e1071)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

source("/home/ijang/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel1.R")


###################################################
#### Load Sanger Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "210937"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "266141"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_exprLayer <- "210931" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "220680" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


# input(continuous), output(<-binay)
classifier<-naiveBayes(iris[,1:4], iris[,5]) 
table(predict(classifier, iris[,-5]), iris[,5])


# random forest


data(iris)
iris <- iris[(iris$Species != "setosa"),]
iris$Species <- factor(iris$Species)
fit <- glm(Species~.,iris,family=binomial)
train.predict <- predict(fit,newdata = iris,type="response")          
library(ROCR)
plot(performance(prediction(train.predict,iris$Species),"tpr","fpr"),col = "red")
auc1 <- performance(prediction(train.predict,iris$Species),"auc")@y.values[[1]]
legend("bottomright",legend=c(paste("Logistic Regression (AUC=",formatC(auc1,digits=4,format="f"),")",sep="")),  
       col=c("red"), lty=1)


library(randomForest)
fit <- randomForest(Species ~ ., data=iris, ntree=50)
train.predict <- predict(fit,iris,type="prob")[,2]          
plot(performance(prediction(train.predict,iris$Species),"tpr","fpr"),col = "red")
auc1 <- performance(prediction(train.predict,iris$Species),"auc")@y.values[[1]]
legend("bottomright",legend=c(paste("Random Forests (AUC=",formatC(auc1,digits=4,format="f"),")",sep="")),  
       col=c("red"), lty=1) 