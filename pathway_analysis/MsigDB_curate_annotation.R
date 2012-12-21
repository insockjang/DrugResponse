# Sample code MsigDB Curation with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
## load the synapse client and login
library(synapseClient)
library(GSA)
library(affy)
synapseLogin()

## set up a project
myName <- "MsigDB Chromosomal Position C1"
ALL <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/msigdb.v3.0.orig.gmt")

C1.POS <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C1/c1.all.v3.0.orig.gmt")

C2.CGP <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C2/Chemical_Genetic_Perturbation/c2.cgp.v3.0.orig.gmt")
C2.BIOCARTA <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C2/Cannonical_Pathway/BIOCARTA/c2.cp.biocarta.v3.0.orig.gmt")
C2.KEGG <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C2/Cannonical_Pathway/KEGG/c2.cp.kegg.v3.0.orig.gmt")
C2.REACTOME <- GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C2/Cannonical_Pathway/REACTOME/c2.cp.reactome.v3.0.orig.gmt")

C3.miRT <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C3/miR_target/c3.mir.v3.0.orig.gmt")
C3.TFT <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C3/TF_target/c3.tft.v3.0.orig.gmt")

C4.CGN <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C4/CancerGene_Neighbor/c4.cgn.v3.0.orig.gmt")
C4.CGM <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C4/CancerGene_Module/c4.cm.v3.0.orig.gmt")

C5.GO_BP <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C5/GO_BP/c5.bp.v3.0.orig(1).gmt")
C5.GO_CC <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C5/GO_CC/c5.cc.v3.0.orig.gmt")
C5.GO_MF <-GSA.read.gmt("/Users/insockjang/Documents/Sage_Challenge/MSigDB/C5/GO_MF/c5.mf.v3.0.orig.gmt")

C1 <- list(POS = C1.POS$geneset.names)
C2 <- list(CGP = C2.CGP$geneset.names, BIOCARTA = C2.BIOCARTA$geneset.names, KEGG = C2.KEGG$geneset.names, REACTOME = C2.REACTOME$geneset.names)
C3 <- list(miRT = C3.miRT$geneset.names, TFT = C3.TFT$geneset.names)
C4 <- list(CGM = C4.CGM$geneset.names, CGN = C4.CGN$geneset.names)
C5 <- list(GO_BP = C5.GO_BP$geneset.names, GO_CC = C5.GO_CC$geneset.names, GO_MF = C5.GO_MF$geneset.names)

save(C1,C2,C3,C4,C5,file = "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_annotation.ROBJECT")

AnnotLayer <- Layer(list(name = "MsigDB annotation", type = "E", parentId = "105191", status="db"))
testFiles <- "/Users/insockjang/Documents/Sage_Challenge/MSigDB/MsigDB_annotation.ROBJECT"
testLayer <- addFile(AnnotLayer, testFiles)
testLayer <- storeEntity(testLayer)
