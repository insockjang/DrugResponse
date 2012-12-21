# Andrew Trister for clinical covariate data
# Sage Bionetworks
# 20120208
#
# This will be used to adjust the clincal covariates from the METABRIC dataset to compare to TCGA for various clinical features

require(Biobase)
require(survival)
require(survcomp)
require(randomSurvivalForest)
require(synapseClient)

synapseLogin()


## PULL IN THE CLINICAL DATA
clinEnt <- loadEntity(120358)
clinDat <- clinEnt$objects$clinical_annotDataFrame
phen <- clinDat@data

#pull in the extended clinical data tables
clinEnt2 <- loadEntity(138645)
clinDat2 <- clinEnt2$objects$clinical_annotDataFrame
phen2 <- clinDat2@data

colnames(phen2)

R1<-rownames(phen)
R2<-rownames(phen2)
phen_01<-data.frame(phen,R1)
phen_02<-data.frame(phen2,R2)

clinAnnotation <- merge(phen_01,phen_02)
clinAnnotation$Time <- clinAnnotation$T/365.25
clinAnnotation <- subset(clinAnnotation, select= -T)
clinAnnotation$T <- ifelse(clinAnnotation$size>20,ifelse(clinAnnotation$size>50,"T3/4","T2" ),"T1")
clinAnnotation$N <- ifelse(clinAnnotation$lymph_nodes_positive>0,1,0)
clinAnnotation$pr_status <- ifelse(clinAnnotation$PR.Expr=="+", "pos", "neg")
clinAnnotation$HER2gain <- ifelse(clinAnnotation$HER2_SNP6_state=="GAIN", "yes", "no")
clinAnnotation$Age35 <- ifelse(clinAnnotation$age_at_diagnosis<35, 1, 0)
clinAnnotation$ERPR <- ifelse(clinAnnotation$er_status=="pos" | clinAnnotation$pr_status=="pos","pos","neg")
clinAnnotation$Histology <- ifelse(clinAnnotation$histological_type=="IDC","Infiltrating Ductal", 
                                   ifelse(clinAnnotation$histological_type=="ILC","Infilitrating Lobular",
                                   ifelse(clinAnnotation$histological_type=="IDC-MED","Medullary",
                                   ifelse(clinAnnotation$histological_type=="IDC-MUC","Mucinous","Mixed Histology"))))
clinAnnotation <- clinAnnotation[!is.na(clinAnnotation$Death),]


#select the features that are in the TCGA set
clinAnnotation <- subset(clinAnnotation, select = c(Histology, age_at_diagnosis, Age35, grade, T, N, ERPR, er_status, pr_status, HER2gain, Time, Death, DeathBreast, Site,R1))

METABRIC_EXP <- clinAnnotation

############################################
# In Sock Jang
# Sage Bionetworks
# MF and CC preprocessing 
# predictive model with cox model
# 2/9/2012

# From clinical covariate data, we only use partial covariates such as Age35, T, N, ERPR, HER2gain, Histology and grade
METABRIC <- METABRIC_EXP[,c("Age35", "T", "N", "ERPR", "HER2gain", "Histology","grade","R1","Time","Death")]
survObj <- Surv(time = METABRIC_EXP[,"Time"], event = METABRIC_EXP[,"Death"])


METABRIC$T[which(METABRIC$T == "T1")] <- 1
METABRIC$T[which(METABRIC$T == "T2")] <- 2
METABRIC$T[-1*union(which(METABRIC$T == "T1"),which(METABRIC$T == "T2"))] <- 3

METABRIC$T<-as.numeric(METABRIC$T)

METABRIC$HER2gain <- ifelse(METABRIC$HER2gain=="yes", 1, 0)

METABRIC$ERPR <- ifelse(METABRIC$ERPR == "pos",1,0)

METABRIC$Histology[which(METABRIC$Histology=="Infilitrating Lobular")] <-1
METABRIC$Histology[which(METABRIC$Histology=="Infiltrating Ductal")] <-2
METABRIC$Histology[which(METABRIC$Histology=="Medullary")] <-3
METABRIC$Histology[which(METABRIC$Histology=="Mixed Histology")] <-4
METABRIC$Histology[which(METABRIC$Histology=="Mucinous")] <-5 
METABRIC$Histology<-as.numeric(METABRIC$Histology)

MMM<-as.matrix(METABRIC[,-8])
rownames(MMM)<-METABRIC$R1
colnames(MMM)<-colnames(METABRIC)[-8]

# Molecular Feature
idExpressionLayer <- "120343"
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

clinicalData <- MMM

dataSets_expr_clinical <- createFeatureAndResponseDataList(exprs(exprData), t(clinicalData))


FEA <- dataSets_expr_clinical$featureData
RES <- dataSets_expr_clinical$responseData


# you can skip above part if you load /home/ijang/Norway_Clinical/feature_clinical.Rdata
# filter NA out from clinical datasets
a<-1:dim(RES)[2]
for(i in 1:dim(RES)[1]){
  a<-intersect(a,which(is.na(RES[i,])==0))
}
RES <- RES[,a]
FEA <- FEA[,a]

training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))

FEA_training <- FEA[,training]
FEA_testing <- FEA[,-training]
RES_training <- RES[,training]
RES_testing <- RES[,-training]

survObj<-Surv(time=RES_training[8,],event = RES_training[9,])

#### modified cox model $1
#### CC + MF
fit_cox1 <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:7),])),survObj,family = "cox")

predTest_cox1 <- predict(fit_cox1,t(rbind(FEA_testing,RES_testing[c(1:7),])),s="lambda.min")
predTrain_cox1 <- predict(fit_cox1,t(rbind(FEA_training,RES_training[c(1:7),])),s="lambda.min")

#### CC only : modified cox model $2 
# This approach means that each clinical covariate has its supplementary molecular features (selected by "feature selection by glmnet")
# Then I combine all features and run cox model
fit_cox2 <-cv.glmnet(t(RES_training[c(1:7),]),survObj,family = "cox")
predTest_cox2 <- predict(fit_cox2,t(RES_testing[c(1:7),]),s="lambda.min")
predTrain_cox2 <- predict(fit_cox2,t(RES_training[c(1:7),]),s="lambda.min")

#### MF only : modified cox model $3
fit_cox3 <-cv.glmnet(t(FEA_training),survObj,family = "cox")
predTest_cox3 <- predict(fit_cox3,t(FEA_testing),s="lambda.min")
predTrain_cox3 <- predict(fit_cox3,t(FEA_training),s="lambda.min")





###### lasso set alpha =1
#### modified cox model $1
#### CC + MF
fit_cox1 <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:7),])),survObj,family = "cox",alpha =1)

predTest_cox1 <- predict(fit_cox1,t(rbind(FEA_testing,RES_testing[c(1:7),])),s="lambda.min")
predTrain_cox1 <- predict(fit_cox1,t(rbind(FEA_training,RES_training[c(1:7),])),s="lambda.min")

#### CC only : modified cox model $2 
# This approach means that each clinical covariate has its supplementary molecular features (selected by "feature selection by glmnet")
# Then I combine all features and run cox model
fit_cox2 <-cv.glmnet(t(RES_training[c(1:7),]),survObj,family = "cox",alpha =1)
predTest_cox2 <- predict(fit_cox2,t(RES_testing[c(1:7),]),s="lambda.min")
predTrain_cox2 <- predict(fit_cox2,t(RES_training[c(1:7),]),s="lambda.min")

#### MF only : modified cox model $3
fit_cox3 <-cv.glmnet(t(FEA_training),survObj,family = "cox",alpha =1)
predTest_cox3 <- predict(fit_cox3,t(FEA_testing),s="lambda.min")
predTrain_cox3 <- predict(fit_cox3,t(FEA_training),s="lambda.min")

library(survcomp)

CONCORD <- function(X, Y, Z) { 
  tt <- concordance.index(x=X, surv.time=Y, surv.event=Z, method="noether", na.rm=TRUE); 
  return(   c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower,"upper"=tt$upper)) 
}


CONCORD(predTrain_cox1,RES_training[10,],RES_training[12,])
CONCORD(predTest_cox1,RES_testing[10,],RES_testing[12,])


CONCORD(predTrain_cox2,RES_training[10,],RES_training[12,])
CONCORD(predTest_cox2,RES_testing[10,],RES_testing[12,])


CONCORD(predTrain_cox3,RES_training[10,],RES_training[12,])
CONCORD(predTest_cox3,RES_testing[10,],RES_testing[12,])


c1<-concordance.index(x=predTrain_cox1,RES_training[10,],RES_training[12,],method ="noether")
c2<-concordance.index(x=predTrain_cox2,RES_training[10,],RES_training[12,],method ="noether")
c3<-concordance.index(x=predTrain_cox3,RES_training[10,],RES_training[12,],method ="noether")
cindex.comp(c1,c2)
cindex.comp(c1,c3)
cindex.comp(c2,c3)

tc1<-concordance.index(x=predTest_cox1,RES_testing[10,],RES_testing[12,],method ="noether")
tc2<-concordance.index(x=predTest_cox2,RES_testing[10,],RES_testing[12,],method ="noether")
tc3<-concordance.index(x=predTest_cox3,RES_testing[10,],RES_testing[12,],method ="noether")
cindex.comp(tc1,tc2)
cindex.comp(tc1,tc3)
cindex.comp(tc2,tc3)
