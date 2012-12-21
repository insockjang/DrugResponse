library(predictiveModeling)
library(synapseClient)
library(survival)
library(survcomp)
library(risksetROC)
library(glmnet)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@") ### not required if configured for automatic login

idExpressionLayer <- "139167"
expressionLayer <- loadEntity(idExpressionLayer)
exprData <- expressionLayer$objects[[1]]

idCopyLayer <- "139169"
copyLayer <- loadEntity(idCopyLayer)
copyData <- copyLayer$objects[[1]]

idClinicalLayer <- "139171"
clinicalLayer <- loadEntity(idClinicalLayer)
clinicalData <- clinicalLayer$objects[[1]]
clinicalData <- clinicalData@data
clinicalData <- clinicalData[!is.na(clinicalData$survDeath), ]
#### prepare feature data for predictive modeling by transposing the matrix to have samples on the rows and features on the columns and scaling the columns
featureData <- scale(t(createAggregateFeatureDataSet(list(expr = exprData, copy = copyData))))

featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "columns")

dataSets_expr_clinical <- createFeatureAndResponseDataList(featureData_filtered, clinicalData)

FEA <- dataSets_expr_clinical$featureData
fea<-t(unique(t(FEA))) # CNV data have lots of replicates

RES <- dataSets_expr_clinical$responseData


# clinical covariate numeric mapping
res<-matrix(0,ncol = 14,nrow = dim(RES)[1])
rownames(res)<-rownames(RES)
colnames(res)<-colnames(RES)[1:14]

res[,"Site"] <- as.numeric(RES[,"Site"])
res[,"ageDiagnosis"]<- as.numeric(RES[,"ageDiagnosis"])
res[,"lymphnodes"] <- ifelse(RES[,"lymphnodes"] == "pos",1,0)
res[,"grade"] <- as.numeric(RES[,"grade"])
res[,"tumorSizeCM"] <- as.numeric(RES[,"tumorSizeCM"])
res[,"chemo"] <- as.numeric(RES[,"chemo"])
res[,"hormone"] <- as.numeric(RES[,"hormone"])
res[,"radiation"] <- as.numeric(RES[,"radiation"])
res[,"HER2"] <- ifelse(RES[,"HER2"]=="pos", 1,0)
res[,"ER"] <- ifelse(RES[,"ER"]=="pos", 1,0)
res[,"PR"] <- ifelse(RES[,"PR"]=="pos", 1,0)
res[,"ERPR"] <- ifelse(RES[,"ERPR"]=="pos", 1,0)
res[,"tripleNegative"] <- ifelse(RES[,"tripleNegative"]=="pos", 1,0)
res[,"histology"] <- ifelse(RES[,"histology"]== "InfilitratingLobular",1,
                            ifelse(RES[,"histology"]== "InfiltratingDuctal",2,
                                   ifelse(RES[,"histology"]== "Medullary",3,
                                          ifelse(RES[,"histology"]== "MixedHistology",4,5))))

rm(RES)
RES<-res

survObj <- Surv(dataSets_expr_clinical$responseData[,"survYears"], dataSets_expr_clinical$responseData[,"survDeath"])

source("/home/ijang/COMPBIO/trunk/users/jang/R5/crossValidationParameterOptimizer_penalty.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/ENetModel.R")

# make Parameter Grids
alpha = seq(0.05,1,0.05)
lambda = createENetTuneGrid(alphas=1)[,2]

foldIndices = createFolds(fea[,1],k=5)

plot_colors <- c("black","blue","red","forestgreen","cyan","magenta")


ALPHA <- foreach(i = 1:length(alpha)) %dopar%{
  aa<-c()
  bb<-c()
  for(k1 in 1:5){
    
    A<-glmnet(fea[-foldIndices[[k1]],],survObj[-foldIndices[[k1]]],alpha= alpha[i],lambda= lambda,family = "cox")
    B<-predict(A,fea[-foldIndices[[k1]],])
    C<-predict(A,fea[foldIndices[[k1]],])
    a<-c()
    b<-c()
    for(k2 in 1:ncol(B)){
      a<-c(a,survConcordance(survObj[-foldIndices[[k1]]] ~ B[,k2])$concordance)
      b<-c(b,survConcordance(survObj[foldIndices[[k1]]] ~ C[,k2])$concordance)      
    }
    
    aa<-rbind(aa,a)
    bb<-rbind(bb,b)     
  }
  
  return(list(Ctrain =aa,Ctest =bb))
}

for(i in 1:20){
  for(k1 in 1:5){
    if(k1 ==1){
      plot(log(sort(lambda,decreasing = TRUE)),ALPHA[[i]]$Ctrain[k1,],ylab = "concordance index",xlab = "log(lambda)",type="b",lty =1,col = plot_colors[k1],ylim = c(0.1,1),pch=k1,main = paste("alpha=",alpha[i],seq=""))
      points(log(sort(lambda,decreasing = TRUE)),ALPHA[[i]]$Ctest[k1,],type="b",lty=2,col = plot_colors[k1],pch=k1)
    }
    else{
      points(log(sort(lambda,decreasing = TRUE)),ALPHA[[i]]$Ctrain[k1,],type="b",col = plot_colors[k1],lty =1,pch=k1)
      points(log(sort(lambda,decreasing = TRUE)),ALPHA[[i]]$Ctest[k1,],type="b",col = plot_colors[k1],lty=2,pch=k1)
    }
  }
}

IQM <- function(vec,na.rm = TRUE){
  a<-quantile(vec,na.rm = TRUE)
  b<-intersect(which(vec>=a[2]),which(vec<=a[4]))
  return(mean(vec[b],na.rm = TRUE))
}

par(mfrow=c(2,1))
plot_colors <- c("black","blue","red","forestgreen","cyan","magenta")

for(i in 1:20){
  if(i ==1){
    plot(log(sort(lambda,decreasing = TRUE)),apply(ALPHA[[i]]$Ctrain,2,IQM),type = "b",pch = i%%6,col=plot_colors[i%%6],lty = i%%6,ylim = c(.4, 1),xlab="log(lambda)",ylab = "cindex",main = "training")
  }
  else{
    points(log(sort(lambda,decreasing = TRUE)),apply(ALPHA[[i]]$Ctrain,2,IQM),type = "b",pch = i%%6,col=plot_colors[i%%6],lty = i%%6)
  }
}

for(i in 1:20){
  if(i ==1){
    plot(log(sort(lambda,decreasing = TRUE)),apply(ALPHA[[i]]$Ctest,2,IQM),type = "b",pch = i%%6,col=plot_colors[i%%6],lty = i%%6,ylim = c(.45, 0.65),xlab="log(lambda)",ylab = "cindex",main = "testing")
  }
  else{
    points(log(sort(lambda,decreasing = TRUE)),apply(ALPHA[[i]]$Ctest,2,IQM),type = "b",pch = i%%6,col=plot_colors[i%%6],lty = i%%6)
  }
}

legend(1,as.character(alpha))
######################################################################################################################################################
#######  Model 2. penalized CC only Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
Model2 <- crossValidationParameterOptimizer_penalty$new()
Model2$findOpt(model=ENetModel, RES, survObj, numFolds = 10,  alpha = alpha, lambda = lambda, family = "cox")

######################################################################################################################################################
#######  Model 3. penalized MF only Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
Model3 <- crossValidationParameterOptimizer_penalty$new()
Model3$findOpt(model=ENetModel, FEA, survObj, numFolds = 10,  alpha = alpha, lambda = lambda, family = "cox")

######################################################################################################################################################
#######  Model 4. penalized CC+MF Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
Model4 <- crossValidationParameterOptimizer_penalty$new()
Model4$findOpt(model=ENetModel, cbind(FEA,RES), survObj, numFolds = 10,  alpha = alpha, lambda = lambda, family = "cox")

######################################################################################################################################################
#######  Model 5. unpenalized CC + penalized MF Model : ENetModel with Cox #######################################################################################
######################################################################################################################################################
Model5 <- crossValidationParameterOptimizer_penalty$new()
Model5$findOpt(model=ENetModel, cbind(FEA,RES), survObj, numFolds = 10,  alpha = alpha, lambda = lambda, family = "cox", penalty.factor = c(rep(1,ncol(FEA)),rep(0,ncol(RES))))


######################################################################################################################################################
#######  submit the Model in order to be evaluated by Sage(Adam)  #######################################################################################
######################################################################################################################################################
submittedModelDatasetId <- "140861"
submittedModelLayer <- Layer(list(name = "In Sock penalized CC only ENet Model", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- addObject(entity=submittedModelLayer, object = Model2)
submittedModelLayer <- storeEntity(submittedModelLayer)

submittedModelDatasetId <- "140861"
submittedModelLayer <- Layer(list(name = "In Sock penalized MF only ENet Model", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- addObject(entity=submittedModelLayer, object = Model3)
submittedModelLayer <- storeEntity(submittedModelLayer)

submittedModelDatasetId <- "140861"
submittedModelLayer <- Layer(list(name = "In Sock penalized CC+MF ENet Model", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- addObject(entity=submittedModelLayer, object = Model4)
submittedModelLayer <- storeEntity(submittedModelLayer)

submittedModelDatasetId <- "140861"
submittedModelLayer <- Layer(list(name = "In Sock unpenalized CC + penalized MF ENet Model", type = "E", parentId = submittedModelDatasetId))
submittedModelLayer <- addObject(entity=submittedModelLayer, object = Model5)
submittedModelLayer <- storeEntity(submittedModelLayer)
# persp function should be used with opt.metricGrid and find optimal alpha and lambda
