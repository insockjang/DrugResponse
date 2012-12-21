# This function is needed in group-lasso-related predictive model 
createAggregateFeatureDataSetGroupIndexForm <- function (dataSetList) 
{
  if (is(dataSetList[[1]], "ExpressionSet")) {
    dataSetList <- lapply(dataSetList, exprs)
  }
  rowNamesList <- lapply(dataSetList, rownames)
  rowAllNames <- Reduce("union",rowNamesList)
      
  colNamesList <- lapply(dataSetList, colnames)
  commonColNames <- Reduce("intersect", colNamesList)
  M<-c()
  for(k in 1:length(dataSetList)){
    m<-match(rowNamesList[[k]],rowAllNames)
    rownames(dataSetList[[k]]) <- paste(m,rowNamesList[[k]], sep = "_")
    M<-c(M,m)
  }
  for (curDsName in names(dataSetList)) {    
    rownames(dataSetList[[curDsName]]) <- sapply(rownames(dataSetList[[curDsName]]), 
                                                 function(rowName) {
                                                   paste(rowName, curDsName, sep = "_")
                                                 }, USE.NAMES = FALSE)
  }
  ds_common_featList <- lapply(dataSetList, function(x) {
    x[, commonColNames]
  })
  ds_feat <- do.call("rbind", ds_common_featList)
  return(list(data = ds_feat, groupIndex = M))
}
