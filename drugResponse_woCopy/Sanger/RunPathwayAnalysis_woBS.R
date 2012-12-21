require(synapseClient)
require(predictiveModeling)

synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/drugResponse/Analysis/pathwayAnalysisFunctions.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  


pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Lasso2/","Lasso",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Ridge2/","Ridge",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/ENet2/","ENet",130,prior = FALSE)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Lasso2/","Lasso",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Ridge2/","Ridge",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/ENet2/","ENet",130,prior = FALSE)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Lasso2/","Lasso",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/Ridge2/","Ridge",130,prior = FALSE)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/COMPBIO/trunk/users/jang/drugResponse_woCopy/Sanger/ENet2/","ENet",130,prior = FALSE)

