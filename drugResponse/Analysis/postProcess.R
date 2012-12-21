
for(kk in 1:24){
  
  filename = paste("cvDrug_",kk,"_Spls.Rdata",sep="")
  load(filename)
  filename1 = paste("cvDrug_",kk,".Rdata",sep="")
  save(resultsRank,resultsScale,file=filename1)
}
  