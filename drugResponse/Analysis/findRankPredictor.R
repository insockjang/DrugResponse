# find order of predictor
findRankPredictor<-function(predictorName,infile,modelName,prior = FALSE,bootstrap = FALSE){
  
  load(infile)
  
  if(!prior){
    if(is.element(modelName, "Ridge")){
      resultsModel<-resultsRidge        
    }
    if(is.element(modelName, "Lasso")){
      resultsModel<-resultsLasso
    }
    if(is.element(modelName, "ENet")){
      resultsModel<-resultsENet
    }
  }
  else{
    if(is.element(modelName, "CGC")){
      resultsModel<-resultsCGC        
    }
    if(is.element(modelName, "TF")){
      resultsModel<-resultsTF
    }
    if(is.element(modelName, "M")){
      resultsModel<-resultsM
    }
    if(is.element(modelName, "R")){
      resultsModel<-resultsR
    }    
  }
  
  if(!bootstrap){
    coefficientsModel<-resultsModel[-1]    
    coefficientModelName<-rownames(resultsModel)[-1]    
    X<-sort(abs(coefficientsModel),decreasing = TRUE,index.return = TRUE)    
    
    sigCoefName<-coefficientModelName[X$ix]
    rankMat<-matrix(NA,ncol=4,nrow=1)
    colnames(rankMat)<-c("expr","copy","mut","all")
    rownames(rankMat)<-predictorName
    
    a<-grep(predictorName,sigCoefName)
    if(length(a)>0){
      for(k in 1:length(a)){
        posterix<-strsplit(sigCoefName[a[k]],"_")[[1]][2]
        prefix<-strsplit(sigCoefName[a[k]],"_")[[1]][1]
        if(!is.na(match(prefix,predictorName))){
          
          if(is.element(posterix,"expr")){
            if(is.na(rankMat[1])){
              rankMat[1]<-a[k]
            }
            else{
              if(rankMat[1]>=a[k]){
                rankMat[1]<-a[k]
              }
            }
          }
          if(is.element(posterix,"copy")){
            if(is.na(rankMat[2])){
              rankMat[2]<-a[k]
            }
            else{
              if(rankMat[2]>=a[k]){
                rankMat[2]<-a[k]
              }
            }
          }
          if(is.element(posterix,"mut")){
            if(is.na(rankMat[3])){
              rankMat[3]<-a[k]
            }
            else{
              if(rankMat[3]>=a[k]){
                rankMat[3]<-a[k]
              }
            }
          }
        }
      }            
    }
    rankMat[4]<-length(sigCoefName)
  }
  # this is for the bootstrapping
  else{
    ResultBS<-c()
    for(k in 1:length(resultsModel)){    
      ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsModel[[k]][-1])),ties.method="min")/length(resultsModel[[k]]))
    }
    rownames(ResultBS)<-rownames(resultsModel[[k]])[-1]  
    reference <- apply(ResultBS,1,sum)        
    X<-sort(reference, decreasing =T, index.return =T)
    sigCoefName<-names(reference)[X$ix]
    
    rankMat<-matrix(NA,ncol=4,nrow=1)
    colnames(rankMat)<-c("expr","copy","mut","all")
    rownames(rankMat)<-predictorName
    
    a<-grep(predictorName,sigCoefName)
    if(length(a)>0){
      for(k in 1:length(a)){
        posterix<-strsplit(sigCoefName[a[k]],"_")[[1]][2]
        prefix<-strsplit(sigCoefName[a[k]],"_")[[1]][1]
        if(!is.na(match(prefix,predictorName))){
          
          if(is.element(posterix,"expr")){
            if(is.na(rankMat[1])){
              rankMat[1]<-a[k]
            }
            else{
              if(rankMat[1]>=a[k]){
                rankMat[1]<-a[k]
              }
            }
          }
          if(is.element(posterix,"copy")){
            if(is.na(rankMat[2])){
              rankMat[2]<-a[k]
            }
            else{
              if(rankMat[2]>=a[k]){
                rankMat[2]<-a[k]
              }
            }
          }
          if(is.element(posterix,"mut")){
            if(is.na(rankMat[3])){
              rankMat[3]<-a[k]
            }
            else{
              if(rankMat[3]>=a[k]){
                rankMat[3]<-a[k]
              }
            }
          }
        }        
      }
    }
    rankMat[4]<-length(sigCoefName)
  }
  return(rankMat)
}
