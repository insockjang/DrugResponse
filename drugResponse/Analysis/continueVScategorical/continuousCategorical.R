require(ROCR)
require(ggplot2)
continuousCategorical<-function(directory,thresholdMethod,numCompound){
  
  for(kk in 1:numCompound){
    filename = paste(directory,"/cvDrug_",kk,".Rdata",sep ="")
    load(filename)
    nfolds = length(resultsScale)
    tePred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testPredictions}
    teObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testObservations}
    
    allTePred<-do.call("c",tePred)
    allTeObsr<-do.call("c",teObsr)
    
    if(is.element(thresholdMethod,"mean")){
      
      allTeObsr1 <-(allTeObsr < mean(allTeObsr))
      allTeObsr2 <-(allTeObsr>= mean(allTeObsr))
      
      Threshold = mean(allTeObsr)
      
      allTeObsr3 <-allTeObsr1 | allTeObsr2
      postTePred<-allTePred[allTeObsr3]
      postTeObsr<-allTeObsr[allTeObsr3]
      postTeObsrFactor<-factor(as.numeric(postTeObsr>=Threshold))
      
    }
    
    if(is.element(thresholdMethod,"mean_sd")){
      
      allTeObsr1 <-(allTeObsr <= (mean(allTeObsr)-0.5*sd(allTeObsr)))
      allTeObsr2 <-(allTeObsr>= (mean(allTeObsr)+ 0.5*sd(allTeObsr)))
      
      Threshold = mean(allTeObsr)
      
      allTeObsr3 <-allTeObsr1 | allTeObsr2
      postTePred<-allTePred[allTeObsr3]
      postTeObsr<-allTeObsr[allTeObsr3]
      postTeObsrFactor<-factor(as.numeric(postTeObsr>=Threshold))
      
    }
    if(is.element(thresholdMethod,"median")){
      allTeObsr1 <-(allTeObsr < median(allTeObsr))
      allTeObsr2 <-(allTeObsr >= median(allTeObsr))
      
      Threshold = median(allTeObsr)
      
      allTeObsr3 <-allTeObsr1 | allTeObsr2
      postTePred<-allTePred[allTeObsr3]
      postTeObsr<-allTeObsr[allTeObsr3]
      postTeObsrFactor<-factor(as.numeric(postTeObsr>=Threshold))
      
    }
    
    
    if(is.element(thresholdMethod,"median_mad")){
      
      allTeObsr1 <-(allTeObsr <= (median(allTeObsr)-0.5*mad(allTeObsr)))
      allTeObsr2 <-(allTeObsr>= (median(allTeObsr)+ 0.5*mad(allTeObsr)))
      
      Threshold = median(allTeObsr)
      
      allTeObsr3 <-allTeObsr1 | allTeObsr2
      postTePred<-allTePred[allTeObsr3]
      postTeObsr<-allTeObsr[allTeObsr3]
      postTeObsrFactor<-factor(as.numeric(postTeObsr>=Threshold))
      
    }
    
    # EVALUATE VALIDATION MODEL PERFORMANCE
    # erPred <- prediction(allTePred,as.numeric(factor(allTeObsr))-1) # factor become c(1,2) from c(0,1) after concatenate
    erPred <- prediction(postTePred,postTeObsrFactor) # factor become c(1,2) from c(0,1) after concatenate
    erPerf <- performance(erPred, "tpr", "fpr")
    erAUC <- performance(erPred, "auc")
    
    
    # FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
    erRFPerf <- performance(erPred, "sens", "spec")
    youdensJ <- erRFPerf@x.values[[1]] + erRFPerf@y.values[[1]] - 1
    jMax <- which.max(youdensJ)
    optCut <- erPerf@alpha.values[[1]][jMax]
    
    optSens <- unlist(erRFPerf@x.values)[jMax]
    optSpec <- unlist(erRFPerf@y.values)[jMax]
    
    #     rankSum <- wilcox.test(validScoreHat[validScore == 0],validScoreHat[validScore == 1])
    
    dfPerf <- as.data.frame(cbind(unlist(erPerf@x.values), unlist(erPerf@y.values)))
    colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")
    
    rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
      geom_line() + 
      geom_abline(slope = 1, colour = "red") +
      opts(title = "Cross Validation ROC Curve") +
      xlab("False Positive Rate") +
      ylab("True Positive Rate") +
      opts(plot.title = theme_text(size = 14))
    
    resultsCat<-list("erPred"=erPred,
                     "erPerf"= erPerf,
                     "erAUC" = erAUC,
                     "rocCurve" = rocCurve)
    
    savename <-paste(directory,"/ROCR_",thresholdMethod,"_",kk,".Rdata",sep ="")
    save(resultsCat,file=savename)
  }
}
