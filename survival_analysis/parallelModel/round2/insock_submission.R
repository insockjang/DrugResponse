###################################################
### loadLibraries
###################################################
library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)


###################################################
### loadData
###################################################
# synapseLogin() ### not required if configured for automatic login

idExpressionLayer <- "160776" ##"160644" 
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idCopyLayer <- "160778" ##"160646"
copyLayer <- loadEntity(idCopyLayer)
copyData <- copyLayer$objects[[1]]

idClinicalFeaturesLayer <- "160780" ##"139171"
clinicalFeaturesLayer <- loadEntity(idClinicalFeaturesLayer)
clinicalFeaturesData <- clinicalFeaturesLayer$objects[[1]]

idClinicalSurvLayer <- "160782"
clinicalSurvLayer <- loadEntity(idClinicalSurvLayer)
clinicalSurvData <- clinicalSurvLayer$objects[[1]]

###################################################
### Call my predictive model class
###################################################

bestCCModel <- setRefClass(Class    = "bestCCModel",                          
                           fields   = c("model"),
                           methods  = list(
                             
                             initialize = function(){                                                                   
                               return(.self)
                             },
                             
                             copy = function() {
                               result <- bestCCModel$new()  				
                               result$model  <- .self$model
                               return(result)
                             },
                             
                             customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...)
                             {                                     
                               RES<-pData(clinicalFeaturesData)
                               # clinical covariate numeric mapping
                               res<-matrix(0,ncol = ncol(RES),nrow = nrow(RES))
                               rownames(res)<-rownames(RES)
                               colnames(res)<-colnames(RES)[1:14]
                               
                               res[,"Site"] <- as.numeric(RES[,"Site"])
                               res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
                               res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
                               res[,"grade"] <- as.numeric(RES[,"grade"])
                               res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                                                           ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                                                  ifelse(RES[,"histology"]== "Medullary",3,
                                                                         ifelse(RES[,"histology"]== "MixedHistology",4,5))))
                               
                               res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
                               res[,"chemo"] <- as.numeric(RES[,"chemo"])
                               res[,"hormone"] <- as.numeric(RES[,"hormone"])
                               res[,"radiation"] <- as.numeric(RES[,"radiation"])
                               res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
                               res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
                               res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
                               res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
                               res[,"tripleNegative"] <- as.numeric(RES[,"tripleNegative"])
                               
                               RES1<-as.data.frame(res)
                               
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               .self$model <- coxph(clinicalSurvData ~ .,RES1[,which(bestCC ==1)])
                               
                               
                             },
                             
                             customPredict = function(exprData, copyData, clinicalFeaturesData){
                               RES<-pData(clinicalFeaturesData)
                               # clinical covariate numeric mapping
                               res<-matrix(0,ncol = ncol(RES),nrow = nrow(RES))
                               rownames(res)<-rownames(RES)
                               colnames(res)<-colnames(RES)[1:14]
                               
                               res[,"Site"] <- as.numeric(RES[,"Site"])
                               res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
                               res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
                               res[,"grade"] <- as.numeric(RES[,"grade"])
                               res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                                                           ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                                                  ifelse(RES[,"histology"]== "Medullary",3,
                                                                         ifelse(RES[,"histology"]== "MixedHistology",4,5))))
                               
                               res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
                               res[,"chemo"] <- as.numeric(RES[,"chemo"])
                               res[,"hormone"] <- as.numeric(RES[,"hormone"])
                               res[,"radiation"] <- as.numeric(RES[,"radiation"])
                               res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
                               res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
                               res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
                               res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
                               res[,"tripleNegative"] <- as.numeric(RES[,"tripleNegative"])
                               
                               RES1<-as.data.frame(res)
                               
                               bestCC <- c(0,1,1,1,0,1,1,1,1,1,0,1,0,0)
                               
                               predict(.self$model, RES1[,which(bestCC ==1)])
                             }
                             )
                           )

###################################################
### trainModel
###################################################
myBestCCModel <- bestCCModel$new()
myBestCCModel$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
trainPredictions <- myBestCCModel$customPredict(exprData, copyData, clinicalFeaturesData)


###################################################
### computeTrainCIndex
###################################################
cIndex_train <- concordance.index(x=trainPredictions, surv.time=clinicalSurvData[,"time"], surv.event=clinicalSurvData[,"status"], na.rm=TRUE, alpha= .05)
cIndex_train$c.index # returns the cindex1
cIndex_train$lower # lower CI bound of cindex1
cIndex_train$upper # upper CI bound of cindex1



###################################################
### submitModel
###################################################
submittedModelParentId <- "161036"
myModelName <- "In Sock best CC model"
submittedModelEntity <- Layer(list(name=myModelName, type="E", parentId=submittedModelParentId))
submittedModelEntity <- addObject(submittedModelEntity, myBestCCModel)
submittedModelEntity <- storeEntity(submittedModelEntity)
