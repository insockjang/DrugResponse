### Pathway analysis with Lasso w/ 100 bootstrapping
library(predictiveModeling)
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("/home/ijang/COMPBIO/trunk/users/jang/pathway_analysis/preRankedTest.R")
source("/home/ijang/COMPBIO/trunk/users/jang/R5/myEnetModel.R")
###################################################
#### Load Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "48339"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$copySet

id_oncomapLayer <- "48341"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$oncomapSet

id_exprLayer <- "48344" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$exprSet

###################################################
### Load Response Data
###################################################

id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, 
                                                  copy = eSet_copy, 
                                                  mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)



for(kk in 1:ncol(dataSets_ccle$responseData)){    
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- t(unique(t(filteredData$featureData)))
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  # set the grids 
  alphas <- c(10^(-5:-1), seq(0.2, 1, by = 0.1))
  lambdas <- exp(seq(-10, 2, length = 30))  
  
  
  # bootstrapping to select the most significant feature set
  coeffLasso <-c()
  boots <- createResample(filteredFeatureDataScaled[,1],times=100)
  for(boot in boots){
    myLassoModel <- myEnetModel$new()    
    
    myLassoModel$customTrain(filteredFeatureDataScaled[boot,], filteredResponseDataScaled[boot], alpha = 1, lambda = lambdas,nfolds = 5)
    
    coeff <- as.vector(myLassoModel$getCoefficients() !=0)
    
    coeffLasso <- cbind(coeffLasso, coeff)
    
    rm(myLassoModel)
  }
  
  # Give manual threshold to select the very significant features.
  threshold = 20
  
  featureNames_all <- colnames(filteredFeatureDataScaled)
  penaltyFactor <- rep(1,length(featureNames_all))
  penaltyFactor[which(apply(coeffLasso,1,sum) >= threshold)] <- 0
  
  finalLassoModel <- myEnetModel$new()      
  finalLassoModel$customTrain(filteredFeatureDataScaled, filteredResponseDataScaled, alpha = 1, lambda = lambdas,nfolds = 5)
  referenceSet <- finalLassoModel$getCoefficients()
  
  
  ref<-referenceSet[-1]
  names(ref) <- rownames(referenceSet)[-1]
  
  a<-sort(abs(ref), decreasing = TRUE, index.return = TRUE)
  
  # select top 200 genes
  referenceGenes<-names(ref)[a$ix[1:200]]
  getUniqueGenesFromFeatureNames <- function(featureNames){
    genes_copy <- sub("_copy","",featureNames[grep("_copy", featureNames)])
    genes_expr <- sub("_expr","",featureNames[grep("_expr", featureNames)])
    genes_mut <- sub("_mut","",featureNames[grep("_mut", featureNames)])  
    geneSet <- union(genes_copy,union(genes_expr, genes_mut))
  }
  
  referenceGenes<-getUniqueGenesFromFeatureNames(referenceGenes)
  
  # MSigDB from synapse
  mSigDB_annotations <- loadEntity(105363)
  mSigDB_symbolID <- loadEntity(105350)
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID  
  allPathways <- mSigDB_annotations$objects$C2$KEGG
  
  # preprocess for FET : make total set
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
  geneAllSet <-c()
  for (i in 1:length(geneAllSetList)){
    geneAllSet<-union(geneAllSet,geneAllSetList[[i]])
  }
  
  analysis <- foreach (i = 1:length(allPathways)) %dopar%{
    print(paste("processing", allPathways[i]))
    
    mSigDB_index <- which(DB$geneset.names == allPathways[i])
    curPathwayGenes <- DB$genesets[[mSigDB_index]]
    
    a1=paste(curPathwayGenes,"_copy",sep ="")
    a2=paste(curPathwayGenes,"_expr",sep ="")
    a3=paste(curPathwayGenes,"_mut",sep ="")
    geneSet <-union(a1,union(a2,a3))
    geneSet <- intersect(geneSet, rownames(referenceSet))
    
    # Here, GSEA part should be taken into consideration
    # to make referenceSet    
    gseaResult <- preRankedTest(abs(r), geneSet,np = 1000)
    
    Mat2x2 <- mat.or.vec(2,2)
    Mat2x2[1,1] <- length(intersect(referenceGenes,geneAllSetList[[i]]))
    Mat2x2[2,1] <- length(setdiff(referenceGenes,geneAllSetList[[i]]))
    Mat2x2[1,2] <- length(setdiff(geneAllSetList[[i]],referenceGenes))
    Mat2x2[2,2] <- length(union(geneAllSet,referenceGenes)) - Mat2x2[1,1]- Mat2x2[1,2]- Mat2x2[2,1]
    fetResult<-fisher.test(Mat2x2)
    return(list(GSEA = gseaResult, FET = fetResult))
  }
  
  filename <- paste("lassoFeatureAnalysisWithBootstrap_",kk,".Rdata",sep ="")
  save(analysis,referenceSet,file=filename)
  
}