# alternative Cox model
# here I tried to predict other clinical informations and with such predicted info, I compared real with predicted survival curve.
# The reason is following.
# With only molecular feature dataset, cox model cannot guaranttee its performance however, other clinical response might be predicted by different models.
# Using these predicted multi responses (such as T, event, er-status, age, etc), we rebuild predicted cox regression model and survival curve might be comparable
# conceptually, molecular feature --> xxx --> hazard ratio
#                             |--------------------|
#
####### Survival
library(predictiveModeling)
library(survival)
library(survcomp)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")

idExpressionLayer <- "120343"
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idClinicalLayer <- "120358"
clinicalLayer <- loadEntity(idClinicalLayer)
clinicalData <- clinicalLayer$objects[[1]]

dataSets_expr_clinical <- createFeatureAndResponseDataList(exprs(exprData), clinicalData)


FEA <- dataSets_expr_clinical$featureData
RES <- dataSets_expr_clinical$responseData

RES1<-matrix(0,nrow = dim(RES)[1],ncol=dim(RES)[2])

RES1[1,]<-as.numeric(RES[1,])
RES1[2,]<-as.numeric(RES[2,])
RES1[3,]<-as.numeric(RES[3,])

RES1[4,which(RES[4,]=="high")]<-1
RES1[4,which(RES[4,]=="moderate")]<-0
RES1[5,]<-as.numeric(RES[5,])

RES1[6,which(RES[6,]=="neg")]<-0
RES1[6,which(RES[6,]=="pos")]<-1

RES1[7,which(RES[7,]=="LOSS")]<-0
RES1[7,which(RES[7,]=="NEUT")]<-1
RES1[7,which(RES[7,]=="GAIN")]<-2
RES1[7,which(RES[7,]=="UNDEF")]<-NA

RES1[8,]<-as.numeric(RES[8,])

RES1[9,which(RES[9,]=="+")]<-1
RES1[9,which(RES[9,]=="-")]<-0

RES1[10,]<-as.numeric(RES[10,])
RES1[11,]<-as.numeric(RES[11,])
RES1[12,]<-as.numeric(RES[12,])

rownames(RES1)<-rownames(RES)
colnames(RES1)<-colnames(RES)

# filter NA out from clinical datasets
a<-1:dim(RES1)[2]
for(i in 1:12){
  a<-intersect(a,which(is.na(RES1[i,])==0))
}
RES1 <- RES1[,a]
FEA <- FEA[,a]

# "T" is transformed to year
Rm<-ceiling(RES1[10,]/30)
Ry<-ceiling(RES1[10,]/365)

training <- sample(1:dim(FEA)[2],round(2*dim(FEA)[2] /3))

FEA_training <- FEA[,training]
FEA_testing <- FEA[,-training]
RES_training <- RES1[,training]
RES_testing <- RES1[,-training]

Rm_training <- Rm[,training]
Rm_testing <- Rm[,-training]

Ry_training <- Ry[,training]
Ry_testing <- Ry[,-training]


# 12 clinical data (responses : ideally this makes cox model fit)
#1 "age_at_diagnosis"        continuous
#2 "grade"                   multinomial 1/2/3
#3 "size"                    continuous 
#4 "cellularity"             binomial high/moderate
#5 "lymph_nodes_positive"    continuous or poisson
#6 "er_status"               binomial
#7 "HER2_SNP6_state"         multinomial GAIN/LOSS/NEUT/UNDEF
#8 "Site"                    multinomial 1/2/3
#9 "PR.Expr"                 binomial +/-
#10 "T"                      continuous(T>0 :day) (it needs to be transformed to month or year) or poisson
#11 "DeathBreast"            binomial 0/1
#12 "Death"                  binomial 0/1

#### modified cox model $1
#### CC + MF
fit_cox1 <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:9),])),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox")

predTest_cox1 <- predict(fit_cox1,t(rbind(FEA_testing,RES_testing[c(1:9),])),s="lambda.min")
predTrain_cox1 <- predict(fit_cox1,t(rbind(FEA_training,RES_training[c(1:9),])),s="lambda.min")

#### CC only : modified cox model $2 
# This approach means that each clinical covariate has its supplementary molecular features (selected by "feature selection by glmnet")
# Then I combine all features and run cox model
fit_cox2 <-cv.glmnet(t(RES_training[c(1:9),]),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox")
predTest_cox2 <- predict(fit_cox2,t(RES_testing[c(1:9),]),s="lambda.min")
predTrain_cox2 <- predict(fit_cox2,t(RES_training[c(1:9),]),s="lambda.min")

#### MF only : modified cox model $3
fit_cox3 <-cv.glmnet(t(FEA_training),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox")
predTest_cox3 <- predict(fit_cox3,t(FEA_testing),s="lambda.min")
predTrain_cox3 <- predict(fit_cox3,t(FEA_training),s="lambda.min")


###### lasso set alpha =1
#### modified cox model $1
#### CC + MF
fit_cox1 <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:9),])),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox",alpha =1)

predTest_cox1 <- predict(fit_cox1,t(rbind(FEA_testing,RES_testing[c(1:9),])),s="lambda.min")
predTrain_cox1 <- predict(fit_cox1,t(rbind(FEA_training,RES_training[c(1:9),])),s="lambda.min")

#### CC only : modified cox model $2 
# This approach means that each clinical covariate has its supplementary molecular features (selected by "feature selection by glmnet")
# Then I combine all features and run cox model
fit_cox2 <-cv.glmnet(t(RES_training[c(1:9),]),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox",alpha =1)
predTest_cox2 <- predict(fit_cox2,t(RES_testing[c(1:9),]),s="lambda.min")
predTrain_cox2 <- predict(fit_cox2,t(RES_training[c(1:9),]),s="lambda.min")

#### MF only : modified cox model $3
fit_cox3 <-cv.glmnet(t(FEA_training),Surv(time=RES_training[10,],event = RES_training[12,]),family = "cox",alpha =1)
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

############################################################################################
############################################################################################
####################   Predict independent Covariate     ###################################
############################################################################################
############################################################################################

fit_cox <-cv.glmnet(t(rbind(FEA_training,RES_training[c(1:9,12),])),Surv(time=RES_training[10,],event = RES_training[11,]),family = "cox")
predTest_cox <- predict(fit_cox,t(rbind(FEA_testing,RES_testing[c(1:9,12),])),s="lambda.min")
predTrain_cox <- predict(fit_cox,t(rbind(FEA_training,RES_training[c(1:9,12),])),s="lambda.min")

CONCORD(predTrain_cox,RES_training[10,],RES_training[11,])
CONCORD(predTest_cox,RES_testing[10,],RES_testing[11,])

a<-createENetTuneGrid()
A<-a[1:2,]
for(i in 1:2){  
  fit <- train(t(FEA_training),RES_training[10,],method = "glmnet",tuneGrid = A[i,],trControl = trainControl(method = "cv"))
}

fit01 <- train(t(FEA_training),RES_training[1,],method = "glmnet")
predTest01 <-predict(fit01,t(FEA_testing))
predTrain01<-predict(fit01,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[1,],predTest01)
plot(RES_training[1,],predTrain01)
dev.off()

fit04 <- train(t(FEA_training),RES_training[4,],method = "glmnet",family = "binomial")
predTest04 <-predict(fit04,t(FEA_testing))
predTrain04<-predict(fit04,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[4,],predTest04)
plot(RES_training[4,],predTrain04)
dev.off()


fit06 <- train(t(FEA_training),RES_training[6,],method = "glmnet",family = "binomial",tuneGrid = createENetTuneGrid(alphas=1))
predTest06 <-predict(fit06,t(FEA_testing))
predTrain06<-predict(fit06,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[6,],predTest06)
plot(RES_training[6,],predTrain06)
dev.off()


fit07 <- train(t(FEA_training),RES_training[7,],method = "glmnet",family = "multinomial",tuneGrid = createENetTuneGrid(alphas=1))
predTest07 <-predict(fit07,t(FEA_testing))
predTrain07<-predict(fit07,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[7,],predTest07)
plot(RES_training[7,],predTrain07)
dev.off()

fit10_p <- train(t(FEA_training),RES_training[10,],method = "glmnet")
predTest10 <-predict(fit10,t(FEA_testing))
predTrain10<-predict(fit10,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[10,],predTest)
plot(RES_training[10,],predTrain)
dev.off()

fit11 <- train(t(FEA_training),RES_training[11,],method = "glmnet",family = "binomial",tuneGrid = createENetTuneGrid(alphas=1))
predTest11 <-predict(fit11,t(FEA_testing))
predTrain11<-predict(fit11,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[11,],predTest11)
plot(RES_training[11,],predTrain11)



library(survival)

binPredTest11 <- predTest11 >=median(predTest11)
binPredTest06 <- predTest06 >=median(predTest06)
binPredTest10 <- predTest10
binPredTest10[which(binPredTest10<=0)] <- 0

survivalObj <- Surv(time = RES1[10,], event= RES1[11,])
survivalObj_realTesting <- Surv(time = RES_testing[10,],event = RES_testing[11,])
survivalObj_predTesting <- Surv(time = binPredTest10, event= binPredTest11)


par(mfrow = c(1,3))
plot(survfit(survivalObj~RES1[6,]),lty= c(1:2),col=c("red","blue"),xlim = c(1,5000))
plot(survfit(survivalObj_realTesting~RES_testing[6,]),lty= c(1:2),col=c("red","blue"),xlim = c(1,5000))
plot(survfit(survivalObj_predTesting~binPredTest06),lty= c(1:2),col=c("red","blue"),xlim = c(1,5000))


# time normalization between 0 and 1
R<-ceiling(RES1[10,]/365)
fit10_year <- train(t(FEA_training),R[training],method = "glmnet",family = "multinomial")
predTest10_year <-predict(fit10_year,t(FEA_testing))
predTrain10_year<-predict(fit10_year,t(FEA_training))

#
y_cox <- Surv(RES_training[10,],RES_training[11,])

fit_cox <- train(t(FEA_training),t(RES_training[10:11,]), method = "glmnet", family = "cox")
predTest_cox <-predict(fit_cox,t(FEA_testing))
predTrain_cox<-predict(fit_cox,t(FEA_training))

##

##
fit10_year_linear <- train(t(FEA_training),Ry[training],method = "glmnet",family = "multinomial",tuneGrid = createENetTuneGrid(alphas=1))
predTest10_year_linear <-predict(fit10_year_linear,t(FEA_testing))
predTrain10_year_linear<-predict(fit10_year_linear,t(FEA_training))

par(mfrow = c(1,2))
plot(RES_testing[11,],predTest10_year_linear)
plot(RES_training[11,],predTrain10_year_linear)
dev.off()

fit10_year_poisson <- train(t(FEA_training),Ry[training],method = "glmnet",family = "poisson",tuneGrid = createENetTuneGrid(alphas=1))
predTest10_year_poisson <-predict(fit10_year_poisson,t(FEA_testing))
predTrain10_year_poisson<-predict(fit10_year_poisson,t(FEA_training))

par(mfrow = c(1,2))
plot(R[-training],predTest10_year_poisson)
plot(R[training],predTrain10_year_poisson)
dev.off()



c1 <- concordance.index(x=RES_testing[6,], surv.time=RES_training[10,], surv.event=RES_training[11,],method="noether")
c2 <- concordance.index(x=P06, surv.time=RES_training[10,], surv.event=RES_training[11,],method="noether")
cindex.comp(c1, c2)


c1 <- concordance.index(x=RES_testing[6,], surv.time=RES_testing[10,], surv.event=RES_testing[11,],method="noether")
c2 <- concordance.index(x=P06, surv.time=RES_testing[10,], surv.event=RES_testing[11,],method="noether")
cindex.comp(c1, c2)