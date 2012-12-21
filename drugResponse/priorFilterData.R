require(predictiveModeling)
priorFilterData <- function(eSet_expr,eSet_copy,eSet_oncomap,adf_drug,priorList = NULL){
  if(!is.null(priorList)){
    
    expr_gene <-rownames(exprs(eSet_expr))
    overlapExprGene<-is.element(expr_gene,priorList)
    
    copy_gene <-rownames(exprs(eSet_copy))
    overlapCopyGene<-is.element(copy_gene,priorList)  
    
    mut_gene <-rownames(exprs(eSet_oncomap))
    overlapMutGene<-is.element(mut_gene,priorList)
    
    eSet_expr<-eSet_expr[overlapExprGene,]
    eSet_copy<-eSet_copy[overlapCopyGene,]
    eSet_oncomap<-eSet_oncomap[overlapMutGene,]
  }
  
  featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))
  
  # NA filter for training set
  featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
  dataSet <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)
  return(dataSet)
}
