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

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/bootstrapLasso2/","Lasso",130)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/bootstrapRidge2/","Ridge",130)

pathwayAnalysisBootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)
pathwayAnalysisBootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/bootstrapENet2/","ENet",130)


pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/Lasso2/","Lasso",130)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/Ridge2/","Ridge",130)

pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"KEGG","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"BIOCARTA","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"REACTOME","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_BP","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_CC","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
pathwayAnalysisNonbootstrap(mSigDB_annotations,DB,"GO_MF","~/PredictiveModel/drugResponse/Sanger/ENet2/","ENet",130)
