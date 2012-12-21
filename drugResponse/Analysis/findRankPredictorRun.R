

library(affy)
source("~/COMPBIO/trunk/users/jang/drugResponse/Analysis/findRankPredictor.R")
id_drugLayer <- "269024" # CCLE drugs 24 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug
drugName_CCLE<-colnames(pData(adf_drug))


id_drugLayer <- "220680" # CCLE drugs 24 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug
drugName_Sanger<-colnames(pData(adf_drug))

drugName <- drugName_CCLE
infile = "~/PredictiveModel/drugResponse_woExp/CCLE/bootstrapLasso2"
# infile = "~/PredictiveModel/drugResponse/CCLE/bootstrapLasso2"
# infileName = paste(infile,"/cvDrug_4.Rdata",sep = "")
# findRankPredictor(predictorName = "BRAF",infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T)

predictor1  ="BRAF"
predictor2  ="NRAS"
predictor3  ="KRAS"
predictor4  ="EGFR"
predictor5  ="ERBB2"
predictor6  ="MDM2"
predictor7  ="TP53"
predictor8  ="RB"
predictor9  ="CDK4"
predictor10 ="CCND1"
predictor11 ="MET"
predictor12 ="ALK"
predictor13 ="FGFR1"
predictor14 ="FGFR2"
predictor15 ="FGFR3"
predictor16 ="IGF1R"
predictor17 ="IRS1"

BRAF<-c()
NRAS<-c()
KRAS<-c()
EGFR<-c()
ERBB2<-c()
MDM2<-c()
TP53<-c()
RB<-c()
CDK4<-c()
CCND1<-c()
MET<-c()
ALK<-c()
FGFR1<-c()
FGFR2<-c()
FGFR3<-c()
IGF1R<-c()
IRS1<-c()

for(k in 1:length(drugName)){
  infileName = paste(infile,"/cvDrug_",k,".Rdata",sep = "")
  BRAF  <-rbind(BRAF,findRankPredictor(predictorName = predictor1,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  NRAS  <-rbind(NRAS,findRankPredictor(predictorName = predictor2,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  KRAS  <-rbind(KRAS,findRankPredictor(predictorName = predictor3,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  EGFR  <-rbind(EGFR,findRankPredictor(predictorName = predictor4,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  ERBB2 <-rbind(ERBB2,findRankPredictor(predictorName = predictor5,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  MDM2  <-rbind(MDM2,findRankPredictor(predictorName = predictor6,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  TP53  <-rbind(TP53,findRankPredictor(predictorName = predictor7,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  RB    <-rbind(RB,findRankPredictor(predictorName = predictor8,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  CDK4  <-rbind(CDK4,findRankPredictor(predictorName = predictor9,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  CCND1 <-rbind(CCND1,findRankPredictor(predictorName = predictor10,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  MET   <-rbind(MET,findRankPredictor(predictorName = predictor11,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  ALK   <-rbind(ALK,findRankPredictor(predictorName = predictor12,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))
  FGFR1 <-rbind(FGFR1,findRankPredictor(predictorName = predictor13,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))  
  FGFR2 <-rbind(FGFR2,findRankPredictor(predictorName = predictor14,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))  
  FGFR3 <-rbind(FGFR3,findRankPredictor(predictorName = predictor15,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))  
  IGF1R <-rbind(IGF1R,findRankPredictor(predictorName = predictor16,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))  
  IRS1 <-rbind(IRS1,findRankPredictor(predictorName = predictor17,infile = infileName, modelName = "Lasso",prior = FALSE, bootstrap = T))  
}
rownames(BRAF)<-drugName
rownames(NRAS)<-drugName
rownames(KRAS)<-drugName
rownames(EGFR)<-drugName
rownames(ERBB2)<-drugName
rownames(MDM2)<-drugName
rownames(TP53)<-drugName
rownames(RB)<-drugName
rownames(CDK4)<-drugName
rownames(CCND1)<-drugName
rownames(MET)<-drugName
rownames(ALK)<-drugName
rownames(FGFR1)<-drugName
rownames(FGFR2)<-drugName
rownames(FGFR3)<-drugName
rownames(IGF1R)<-drugName
rownames(IRS1)<-drugName

save(BRAF,NRAS,KRAS,EGFR,ERBB2,MDM2,TP53,RB,CDK4,CCND1,MET,ALK,FGFR1,FGFR2,FGFR3,IGF1R,IRS1,file=paste(infile,"/predictorRank.Rdata",sep=""))
