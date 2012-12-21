require(synapseClient)
require(predictiveModeling)
#require(multicore)
#require(doMC)
#registerDoMC()
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
source("~/COMPBIO/trunk/users/jang/pathway_analysis/pathwayAnalysisFunctions.R")
################################################################## bs analysis
# MSigDB from synapse
mSigDB_annotations <- loadEntity("syn105363")
mSigDB_symbolID <- loadEntity("syn105350")
DB<-mSigDB_symbolID$objects$MsigDB_symbolID  

# pathwayAnalysisBootstrap(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE)
# pathwayAnalysisBootstrap<-function(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE)
# pathwayAnalysisNonbootstrap<-function(mSigDB_annotations,DB,pathwayName,bootstrapFolder,modelName,numCompounds,prior = FALSE)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2/","Lasso",24)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/bootstrapRidge2/","Ridge",24)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/bootstrapENet2/","ENet",24)


pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/Lasso2/","Lasso",24)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/Ridge2/","Ridge",24)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/CCLE/ENet2/","ENet",24)
