# This function is needed in group-lasso-related predictive model 
groupIndexForm <- function (dataMat) {
  colNames<-sapply(colnames(dataMat), 
                   function(rowName) {
                      a<-strsplit(rowName, "_")
                      return(a[[1]][[1]])
                    }, 
                   USE.NAMES = FALSE)
  
  colAllNames <- unique(colNames)
  
  groupIndex<-match(colNames,colAllNames)
  return(groupIndex)
  }
