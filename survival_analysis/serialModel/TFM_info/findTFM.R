library(affy)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)

library(synapseClient)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
#preprocessed results from lasso in order to select Molecular features

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

findTFM <- function(eSet){
  x <- lumiHumanAllENTREZID
  TF<-read.delim("/Volumes/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TFM_info/tf_list_nov08_entrez_1875.txt",header = F)
  M<-read.table("/Volumes/ijang/COMPBIO/trunk/users/jang/survival_analysis/serialModel/TFM_info/modulator_list.txt",header = F)
  # Get the probe identifiers that are mapped to an ENTREZ Gene ID
  mapped_probes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_probes])
  
  probeGeneMap <-c()
  for(i in 1:length(xx)){
    probeGeneMap <- rbind(probeGeneMap,c(as.character(names(xx[i])),as.integer(xx[[i]])))
  }
  
  featureID <-featureNames(eSet)
  illuminaID <- IlluminaID2nuID(featureID, lib='lumiHumanIDMapping')
  
  a1<-is.element(probeGeneMap[,2],TF[,1])
  a2<-is.element(probeGeneMap[,2],M[,1])
  
  b1<-is.element(illuminaID[,7],probeGeneMap[a1,1]) # TF
  b2<-is.element(illuminaID[,7],probeGeneMap[a2,1]) # M
  return(list(TF=b1,M=b2))
}

exprTFM<-findTFM(exprData)
copyTFM<-findTFM(copyData)
