
fileFolder <- "/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/prior/CCLE/Ridge/"
fileFolder1 <- "/Volumes/ijang/COMPBIO/trunk/users/jang/drugResponse/CCLE/Ridge/"

corPearson<-c()

for(kk in 1:24){
  
  filename <- paste(fileFolder,"cvDrug_",kk,".Rdata",sep="")
  
  load(filename)
  trPred <- foreach(k = 1:5) %do%{resultsCGC[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{resultsCGC[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{resultsCGC[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{resultsCGC[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  a1<-cor(allTePred,allTeObsr)
  
  trPred <- foreach(k = 1:5) %do%{resultsTF[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{resultsTF[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{resultsTF[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{resultsTF[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  a2<-cor(allTePred,allTeObsr)
  
  trPred <- foreach(k = 1:5) %do%{resultsM[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{resultsM[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{resultsM[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{resultsM[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  a3<-cor(allTePred,allTeObsr)
  
  trPred <- foreach(k = 1:5) %do%{resultsR[[k]]$trainPredictions}
  tePred <- foreach(k = 1:5) %do%{resultsR[[k]]$testPredictions}
  trObsr <- foreach(k = 1:5) %do%{resultsR[[k]]$trainObservations}
  teObsr <- foreach(k = 1:5) %do%{resultsR[[k]]$testObservations}
  
  allTrPred<-do.call("c",trPred)
  allTePred<-do.call("c",tePred)
  allTrObsr<-do.call("c",trObsr)
  allTeObsr<-do.call("c",teObsr)
  
  a4<-cor(allTePred,allTeObsr)
  
  corPearson<-rbind(corPearson,c(a1,a2,a3,a4))
  
}

filename1 <- paste(fileFolder,"cv5performance.Rdata",sep="")

save(corPearson,file = filename1)

load(paste(fileFolder1,"cv5performance.Rdata",sep=""))

plot(corPearsonScale[,2])
points(corPearson[,1],col="red")
points(corPearson[,2],col="blue")
points(corPearson[,3],col="green")
points(corPearson[,4],col="cyan")
