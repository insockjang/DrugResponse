crossValidatePredictiveModelLogistic <- 
  function(featureData, responseData, model, numFolds = 5, ...){
    
    #-----------------------------------------------------------------------
    # Split the data into training and test partitions
    # -----------------------------------------------------------------------
    set.seed(2)
    foldIndices <- createFolds(featureData[,1], k = numFolds, list = TRUE)
    
    message("Training partitions")
    
    foldResults <- foreach(fold = foldIndices) %dopar% {
      roc<-function(A){
        a<-matrix(0,nrow=2,ncol=2)
        for(kkk in 1:nrow(A)){
          if(A[kkk,1]==1 && A[kkk,2]==1) a[1,1]=a[1,1]+1
          if(A[kkk,1]==1 && A[kkk,2]==0) a[1,2]=a[1,2]+1
          if(A[kkk,1]==0 && A[kkk,2]==1) a[2,1]=a[2,1]+1
          if(A[kkk,1]==0 && A[kkk,2]==0) a[2,2]=a[2,2]+1    
        }        
        TP<-a[1,1]
        FP<-a[1,2]
        FN<-a[2,1]
        TN<-a[2,2]
        TPR <-TP/(TP+FN)
        FPR <-FP/(FP+TN)  
        b<-c(TPR,FPR)
        names(b)<-c("TPR","FPR")
        return(b)
      }
      trResponse<-responseData[-fold]
      teResponse<-responseData[fold]
      
      trFeature<-featureData[-fold,]
      teFeature<-featureData[fold,]
      
      b<-sort(trResponse)
      threshold<-b[(seq_along(b)+2)%%5==0]
      #print(threshold)
      trFeature1 <-as.data.frame(trFeature)
      teFeature1 <-as.data.frame(teFeature)
      
      #print(min(trResponse))
      ROC<-c()
      for(kk in 2:length(threshold)-1){
        
        trResponse1 <- factor(trResponse<=threshold[kk],levels = c("FALSE","TRUE"))
        model$customTrain(trFeature1, trResponse1, ...)
        pred<-model$customPredict(teFeature)
        obs <- teResponse<=threshold[kk]
        #print(rbind(as.numeric(pred==1),as.numeric(obs==1)))
        ROC<-rbind(ROC,roc(cbind(as.numeric(as.logical.factor(pred)==1),as.numeric(obs==1))))
      }
      print(ROC)
      return(ROC)           
    }
    
    results<-c()
    for(k in 1:length(foldResults)){
      results<-rbind(results,foldResults[[k]])
    }
    return(results)
  }
