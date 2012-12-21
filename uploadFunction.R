
uploadSynapseMsigDB<-function(featureSets,geneSetDbs,algorithms){
  source("~/COMPBIO/trunk/users/jang/myFunctions.R")
  id_drugLayer <- "269024" # CCLE drugs 24 
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  
  curFeatureSet<-featureSets
  if(!is.na(match(curFeatureSet,"exprCopyMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/"
  if(!is.na(match(curFeatureSet,"copyMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woExp/CCLE/"
  if(!is.na(match(curFeatureSet,"exprMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/"
  if(!is.na(match(curFeatureSet,"exprCopy"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woMut/CCLE/"
  
  curAlgorithm<-algorithms
  curGeneSetDb<-geneSetDbs
  if(is.element(curGeneSetDb,"GO_BP")) curPathways <- mSigDB_annotations$objects$C5$GO_BP
  if(is.element(curGeneSetDb,"GO_CC")) curPathways <- mSigDB_annotations$objects$C5$GO_CC
  if(is.element(curGeneSetDb,"GO_MF")) curPathways <- mSigDB_annotations$objects$C5$GO_MF
  if(is.element(curGeneSetDb,"KEGG")) curPathways <- mSigDB_annotations$objects$C2$KEGG
  if(is.element(curGeneSetDb,"REACTOME")) curPathways <- mSigDB_annotations$objects$C2$REACTOME
  if(is.element(curGeneSetDb,"BIOCARTA")) curPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  
  for (drugIndex in 1:length(drugName)){
    
    curDrugName <- drugName[drugIndex]
    
    filePath <- paste(filePath1, curAlgorithm, "/", curGeneSetDb, "_", drugIndex, ".Rdata", sep="")
    print(filePath)
    load(filePath)
    for (pathwayIndex in 1:length(curPathways)){
      
      curPathway <- make.names(curPathways[pathwayIndex])
      curAnalyticResult <- analyticResult[[pathwayIndex]]
      
      entityName <- paste(curDrugName, curPathway, curAlgorithm, curFeatureSet, curGeneSetDb, sep="_")
      
      #       queryEntity <- paste("SELECT id FROM entity WHERE entity.name == \"",entityName,"\"",sep="")
      #       queryResults <- synapseQuery(queryEntity)
      #       
      #       if(is.null(queryResults)){
      newEntity <- Data(list(name= entityName, parentId = "syn582135"))          
      newEntity <- createEntity(newEntity)          
      
      newEntity$properties$name <- entityName
      
      newEntity$annotations$drug <- curDrugName
      newEntity$annotations$pathway <- curPathway
      
      newEntity$annotations$featureSet <- curFeatureSet
      newEntity$annotations$geneSetDb <- curGeneSetDb
      newEntity$annotations$algorithm <- curAlgorithm
      
      if(is.nan(curAnalyticResult$gseaResult$nes)){
        newEntity$annotations$gseaEs = "---"
        newEntity$annotations$gseaNes = "---"
        newEntity$annotations$gseaPValue = "---" 
      } else {
        newEntity$annotations$gseaEs = curAnalyticResult$gseaResult$es
        newEntity$annotations$gseaNes = curAnalyticResult$gseaResult$nes
        newEntity$annotations$gseaPValue = curAnalyticResult$gseaResult$p.value
        if(curAnalyticResult$gseaResult$p.value <= 0.05){
          addObject(newEntity, curAnalyticResult)  
        }
      }
      
      #             newEntity$annotations$fetOddsRatio = curAnalyticResult$fetResult$estimate
      if(is.finite(curAnalyticResult$fetResult$estimate)){
        newEntity$annotations$fetOddsRatio = curAnalyticResult$fetResult$estimate
      }
      else{
        newEntity$annotations$fetOddsRatio = "INF"
      } 
      
      newEntity$annotations$fetPValue = curAnalyticResult$fetResult$p.value
      
      storeEntity(newEntity)
      #       }
    }
  }
}


uploadSynapseGraphite<-function(featureSets,geneSetDbs,algorithms){
  source("~/COMPBIO/trunk/users/jang/myFunctions.R")
  id_drugLayer <- "269024" # CCLE drugs 24 
  layer_drug <- loadEntity(id_drugLayer)
  adf_drug <- layer_drug$objects$adf_drug
  drugName<-colnames(pData(adf_drug))
  
  curFeatureSet<-featureSets
  if(!is.na(match(curFeatureSet,"exprCopyMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse/CCLE/"
  if(!is.na(match(curFeatureSet,"copyMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woExp/CCLE/"
  if(!is.na(match(curFeatureSet,"exprMut"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/"
  if(!is.na(match(curFeatureSet,"exprCopy"))) filePath1 <- "~/COMPBIO/trunk/users/jang/drugResponse_woMut/CCLE/"
  
  curAlgorithm<-algorithms
  curGeneSetDb <- geneSetDbs
  if(is.element(curGeneSetDb,"KEGG")) curPathways <- KEGG
  if(is.element(curGeneSetDb,"BIOCARTA")) curPathways <- BIOCARTA
  if(is.element(curGeneSetDb,"REACTOME")) curPathways <- REACTOME
  if(is.element(curGeneSetDb,"NCI")) curPathways <- NCI
  
  for (drugIndex in 1:length(drugName)){
    curDrugName <- drugName[drugIndex]
    
    filePath <- paste(filePath1, curAlgorithm, "/", curGeneSetDb, "_graphite_", drugIndex, ".Rdata", sep="")
    print(filePath)
    load(filePath)
    for (pathwayIndex in 1:length(curPathways)){
      
      curPathway <- make.names(curPathways[[pathwayIndex]]@title)
      curPathway <- gsub("ß","beta",curPathway)
      
      curAnalyticResult <- analyticResult[[pathwayIndex]]
      
      entityName <- paste(curDrugName, curPathway, curAlgorithm, curFeatureSet, curGeneSetDb, sep="_")
      #       queryEntity <- paste("SELECT id FROM entity WHERE entity.name == \"",entityName,"\"",sep="")
      #       queryResults <- synapseQuery(queryEntity)
      #       if(is.null(queryResults)){
      newEntity <- Data(list(name= entityName, parentId = "syn582134"))          
      newEntity<-createEntity(newEntity)          
      newEntity$annotations$drug <- curDrugName
      newEntity$annotations$pathway <- curPathway          
      newEntity$annotations$featureSet <- curFeatureSet
      newEntity$annotations$geneSetDb <- curGeneSetDb
      newEntity$annotations$algorithm <- curAlgorithm
      
      if(is.nan(curAnalyticResult$gseaResult$nes)){
        newEntity$annotations$gseaEs = "---"
        newEntity$annotations$gseaNes = "---"
        newEntity$annotations$gseaPValue = "---" 
      } else {
        newEntity$annotations$gseaEs = curAnalyticResult$gseaResult$es
        newEntity$annotations$gseaNes = curAnalyticResult$gseaResult$nes
        newEntity$annotations$gseaPValue = curAnalyticResult$gseaResult$p.value  
        if(curAnalyticResult$gseaResult$p.value <= 0.05){
          addObject(newEntity, curAnalyticResult)
        }              
      }
      if(is.finite(curAnalyticResult$fetResult$estimate)){
        newEntity$annotations$fetOddsRatio = curAnalyticResult$fetResult$estimate
      }
      else{
        newEntity$annotations$fetOddsRatio = "INF"
      } 
      newEntity$annotations$fetPValue = curAnalyticResult$fetResult$p.value            
      storeEntity(newEntity)
      }
    #     }
  }
}
