require(synapseClient)
require(predictiveModeling)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/drugResponse/Analysis/pathwayAnalysisFunctions.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  


pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapLasso2/","Lasso",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapRidge2/","Ridge",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapENet2/","ENet",24,prior = FALSE)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapLasso2/","Lasso",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapRidge2/","Ridge",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapENet2/","ENet",24,prior = FALSE)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapLasso2/","Lasso",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapRidge2/","Ridge",24,prior = FALSE)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/CCLE/bootstrapENet2/","ENet",24,prior = FALSE)

