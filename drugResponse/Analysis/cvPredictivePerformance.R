cvPredictivePerformance<-function(dirNames,kk,prior = FALSE){
  
  indTestCor<-c()
  indTrainCor<-c()
  overallTestCor<-c()
  overallTrainCor<-c()
  
  for(kk in 1:kk){
    
    filename = paste(dirNames,"/cvDrug_",kk,".Rdata",sep="")
    if(is.element(prior,"CGC")){
      resultsScale <- resultsCGC      
    }
    if(is.element(prior,"TF")){
      resultsScale <- resultsTF      
    }
    if(is.element(prior,"M")){
      resultsScale <- resultsM      
    }
    if(is.element(prior,"R")){
      resultsScale <- resultsR      
    }
    
    load(filename)
    
    testCor <- c()
    trainCor <- c()
    
    for(k in 1:5){
      testCor <- c(testCor,cor(resultsScale[[k]]$testPredictions,resultsScale[[k]]$testObservations))
      trainCor <- c(testCor,cor(resultsScale[[k]]$trainPredictions,resultsScale[[k]]$trainObservations))
    }
    
    indTestCor<-cbind(indTestCor,testCor)
    indTrainCor<-cbind(indTrainCor,trainCor)
    
    trPred <- foreach(k = 1:5) %do%{resultsScale[[k]]$trainPredictions}
    tePred <- foreach(k = 1:5) %do%{resultsScale[[k]]$testPredictions}
    trObsr <- foreach(k = 1:5) %do%{resultsScale[[k]]$trainObservations}
    teObsr <- foreach(k = 1:5) %do%{resultsScale[[k]]$testObservations}
    
    allTrPred<-do.call("c",trPred)
    allTePred<-do.call("c",tePred)
    allTrObsr<-do.call("c",trObsr)
    allTeObsr<-do.call("c",teObsr)
    
    allTestCor<-cor(allTePred,allTeObsr)
    allTrainCor<-cor(allTrPred,allTrObsr)
    
    overallTestCor<-c(overallTestCor,allTestCor)
    overallTrainCor<-c(overallTrainCor,allTrainCor)
  }
  return(list(indTestCor=indTestCor,indTrainCor=indTrainCor,overallTestCor=overallTestCor,overallTrainCor=overallTrainCor))
}
