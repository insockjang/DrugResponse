require(synapseClient)
require(predictiveModeling)
require(multicore)
require(doMC)
require()
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/drugResponse/Analysis/pathwayAnalysisFunctions.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  

# pathwayAnalysisBootstrap(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE)
  
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse/prior/CCLE/bootstrapLasso2","CGC",24,prior = TRUE)
