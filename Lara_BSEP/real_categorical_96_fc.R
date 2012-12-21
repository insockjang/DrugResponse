### LARA BSEP classifier 
library(predictiveModeling)
library(synapseClient)
library(ggplot2)
library(ROCR)
library(modeest)
require(multicore)
require(doMC)
registerDoMC()
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/R5/crossValidatePredictiveModel_categorical.R")
source("~/COMPBIO/trunk/users/jang/R5/myCatEnetModel.R")
source("~/COMPBIO/trunk/users/jang/R5/myRandomForestModel.R")

###################################################
#### Load BSEP Molecular Feature Data from Synapse ####
###################################################

data_expr <- train96ld
Name<-names(train96catBSEP)
data_drug<-as.numeric(as.matrix(train96catBSEP))
names(data_drug)<-Name


# NA filter for training set
featureData_filtered <- filterNasFromMatrix(data_expr, filterBy = "rows")
dataSets_bsep <- createFeatureAndResponseDataList(t(featureData_filtered),data_drug)

filteredData<-filterPredictiveModelData(dataSets_bsep$featureData,dataSets_bsep$responseData[,drop=FALSE])
filteredFeatureData  <- filteredData$featureData
filteredResponseData <- filteredData$responseData


alphas = seq(0.1,1,by = 0.1)
lambdas = createENetTuneGrid(alphas = 1)[,2]


modelLassoCat<-myCatEnetModel$new()
modelLassoCat$customTrain(filteredFeatureData,filteredResponseData,alpha=1,lambda = lambdas, nfolds=5)


allTePred <- modelLassoCat$customPredict(t(as.matrix(test96ld)))
allTeObsr <- as.numeric(test96catBSEP)

a1<-which(!is.na(allTePred)==1)
a2<-which(!is.na(allTeObsr)==1)
a<-intersect(a1,a2)

AllTePred<-allTePred[a]
AllTeObsr<-allTeObsr[a]

# EVALUATE VALIDATION MODEL PERFORMANCE
erPred <- prediction(AllTePred,AllTeObsr-1) # factor become c(1,2) from c(0,1) after concatenate
erPerf <- performance(erPred, "tpr", "fpr")
erAUC <- performance(erPred, "auc")

# FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
erRFPerf <- performance(erPred, "sens", "spec")
youdensJ <- erRFPerf@x.values[[1]] + erRFPerf@y.values[[1]] - 1
jMax <- which.max(youdensJ)
optCut <- erPerf@alpha.values[[1]][jMax]

optSens <- unlist(erRFPerf@x.values)[jMax]
optSpec <- unlist(erRFPerf@y.values)[jMax]

#     rankSum <- wilcox.test(validScoreHat[validScore == 0],validScoreHat[validScore == 1])

dfPerf <- as.data.frame(cbind(unlist(erPerf@x.values), unlist(erPerf@y.values)))
colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")

rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
  geom_line() + 
  geom_abline(slope = 1, colour = "red") +
  opts(title = "Cross Validation ROC Curve") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  opts(plot.title = theme_text(size = 14))

return(list("erPred"=erPred,
            "erPerf"=erPerf,
            "erAUC"=erAUC,
            "rocCurve" = rocCurve))

resultsLasso<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=1, lambda = lambdas)
resultsENet<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myCatEnetModel$new(), alpha=alphas, lambda = lambdas)
resultsRandomForest_300<-crossValidatePredictiveModel_categorical(filteredFeatureData, filteredResponseData, model = myRandomForestModel$new(), numFolds = 5,ntree = 300)


filename = paste("~/COMPBIO/trunk/users/jang/Lara_BSEP/ROCR_categorical_cv5_96_fc.Rdata",sep = "")
save(resultsLasso,resultsENet,resultsRandomForest_300,file = filename)

