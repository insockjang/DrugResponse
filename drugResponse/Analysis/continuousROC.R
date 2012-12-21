nfolds = 5
for(kk in 1:130){
  filename = paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/Sanger/SVM/cvDrug_",kk,".Rdata",sep ="")
  load(filename)
  trPred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainPredictions}
  tePred <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testPredictions}
  trObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$trainObservations}
  teObsr <- foreach(k = 1:nfolds) %do%{resultsScale[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  Threshold = median(allTeObsr)
  allTeObsr1 <-(allTeObsr <= (median(allTeObsr)-0.5*mad(allTeObsr)))
  allTeObsr2 <-(allTeObsr>= (median(allTeObsr)+ 0.5*mad(allTeObsr)))
  
  allTeObsr3 <-allTeObsr1 | allTeObsr2
  postTePred<-allTePred[allTeObsr3]
  postTeObsr<-allTeObsr[allTeObsr3]
  postTeObsrFactor<-factor(as.numeric(postTeObsr>=Threshold))
  
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
  savename <-paste("/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/Sanger/SVM/ROCR_median_MAD_",kk,".Rdata",sep ="")
  save(resultsCat,file=savename)
}
