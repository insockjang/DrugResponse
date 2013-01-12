library(synapseClient)
library(predictiveModeling)
library(multicore)
source("~/DrugResponse/R5/bootstrapPredictiveModel.R")
source("~/DrugResponse/R5/myEnetModel.R")

id_exprLayer <- "syn1532991" 
layer_expr <- loadEntity(id_exprLayer)
res2 <- layer_expr$objects$res2

id_headerLayer <- "syn1532987" 
layer_header <- loadEntity(id_headerLayer)
header <- layer_header$objects$header

id_proteinLayer <- "syn1583864" # "syn1583864"
layer_protein <- loadEntity(id_proteinLayer)
protein2 <- layer_protein$objects$protein_res2

id_headerLayer <- "syn1583866"
layer_header <- loadEntity(id_headerLayer)
header2 <- layer_header$objects$header_trimmed
colnames(protein2)<-as.character(header2$SAMPLE.ID)

featureData <- createAggregateFeatureDataSet(list(expr = res2, prot = protein2))

data_drug <- as.numeric(as.matrix(header2$BRAAK))
names(data_drug)<-as.character(header2$SAMPLE.ID)


dataSets <- createFeatureAndResponseDataList(t(featureData),data_drug)


TRIX<-as.numeric(as.matrix(header$TRIX))
names(TRIX)<- as.character(header$SAMPLE.ID)
trix<-TRIX[rownames(dataSets$featureData)]


# MCI
aa<-which(trix == 2)
aaa<-which(is.na(dataSets$responseData)==0)                                        
AA<-intersect(aa,aaa)


filteredFeatureData  <- scale(dataSets$featureData[AA,])
filteredResponseData <- scale(dataSets$responseData[AA])

filteredFeatureData_sparse<-filteredFeatureData
filteredFeatureData_sparse[which(is.na(filteredFeatureData)==1)]<-0

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]


resultsLasso<-bootstrapPredictiveModel(filteredFeatureData_sparse, filteredResponseData, model = myEnetModel$new(), numBootstrap=100, alpha=1, lambda = lambdas)
resultsENet<-bootstrapPredictiveModel(filteredFeatureData_sparse, filteredResponseData, model = myEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)

save(resultsLasso,resultsENet,file = "~/Results/AD/BRAAK/MCI/bs100Results_exprProtein2.Rdata")




id_exprLayer <- "syn1532991" 
layer_expr <- loadEntity(id_exprLayer)
res2 <- layer_expr$objects$res2

id_headerLayer <- "syn1532987" 
layer_header <- loadEntity(id_headerLayer)
header <- layer_header$objects$header

id_proteinLayer <- "syn1583862" # "syn1583864"
layer_protein <- loadEntity(id_proteinLayer)
protein2 <- layer_protein$objects$protein_res2

id_headerLayer <- "syn1583866"
layer_header <- loadEntity(id_headerLayer)
header2 <- layer_header$objects$header_trimmed
colnames(protein2)<-as.character(header2$SAMPLE.ID)

featureData <- createAggregateFeatureDataSet(list(expr = res2, prot = protein2))

data_drug <- as.numeric(as.matrix(header2$BRAAK))
names(data_drug)<-as.character(header2$SAMPLE.ID)


dataSets <- createFeatureAndResponseDataList(t(featureData),data_drug)


# MMSE
TRIX<-as.numeric(as.matrix(header$TRIX))
names(TRIX)<- as.character(header$SAMPLE.ID)
trix<-TRIX[rownames(dataSets$featureData)]


# MCI
aa<-which(trix == 2)
aaa<-which(is.na(dataSets$responseData)==0)                                        
AA<-intersect(aa,aaa)


filteredFeatureData  <- scale(dataSets$featureData[AA,])
filteredResponseData <- scale(dataSets$responseData[AA])

filteredFeatureData_sparse<-filteredFeatureData
filteredFeatureData_sparse[which(is.na(filteredFeatureData)==1)]<-0

alphas =unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]


resultsLasso<-bootstrapPredictiveModel(filteredFeatureData_sparse, filteredResponseData, model = myEnetModel$new(), numBootstrap=100, alpha=1, lambda = lambdas)
resultsENet<-bootstrapPredictiveModel(filteredFeatureData_sparse, filteredResponseData, model = myEnetModel$new(), numBootstrap=100, alpha=alphas, lambda = lambdas)

save(resultsLasso,resultsENet,file = "~/Results/AD/BRAAK/MCI/bs100Results_exprProtein1.Rdata")
