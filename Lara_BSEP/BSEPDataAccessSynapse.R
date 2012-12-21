#########################################################################################
##Loading the training and test data for the classification of BSEP inhibition project.##
##Lara Mangravite, November 15 2012                                                    ##
#########################################################################################

#All data is stored in Synapse in the project: Toxicological Classification of Compounds
#which can be viewed at:  https://synapse.sagebase.org/#Synapse:syn375898
#this is a private project: To gain access, you will need to create a Synapse account and then ask me to add you to this project.

##Load Synapse Client.
source('http://depot.sagebase.org/CRAN.R')
pkgInstall(c("synapseClient"))
library(synapseClient)
synapseCacheDir('~/synapseCache/')
#synapseLogin("YOUR_SYNAPSE_LOGIN_ID_HERE") #EG, synapseLogin("lara.mangravite@sagebase.org")
#enter your Synapse password to continue

########################TRANING DATA#####################################################
##Training data is located in folder: syn480538
##Data was transfered from the Japanese Toxicology project (http://toxico.nibio.go.jp/open-tggates/search.html)
##Data was normalized using SNM and Drug Response signatures were calculated separately for compound 
#by calculating either the fold change (H/aveC) or Logged Difference (logH - logaveC)
#for each expression signature following high dose exposure (H) relative the average expression for the three control-exposed samples (aveC)

qryTRAIN <- synapseQuery('select id, name from entity where entity.parentId=="syn480538"')

#####################################
##TRAINING DATA FOR 96-HOUR EXPOSURES
#####################################
#Contains 236 drug response signatures, representing 79 compounds (9 BSEP inhibitors and 65 noninhibitors)

#Expression response to drug exposure
train96<-loadEntity(qryTRAIN[7,2])
train96fc<-train96$objects$FoldChangeResponse  #Drug response signature calculated as fold change (H/aveC)
train96ld<-train96$objects$LogDeltaResponseData #Drug response signature calculated as difference in logs (logH - logaveC)
probeMap<-train96$objects$ProbeAnnotations
#Phenotypes
trainp96<-loadEntity(qryTRAIN[5,2])
train96p<-trainp96$objects$Phenotypes  #BSEP phenotypes
train96m<-trainp96$objects$Metadata  #information about originating CEL files and the rats.

#Continuous trait for BSEP inhibition (IC50 for rat BSEP inhibition using in vitro assay)
train96IC50<-train96p$RatBSEPIC50
names(train96IC50)<-train96p$SampleID
#Categorical trait for BSEP inhibition.  0 where ratIC50>=100 and 1 where ratIC50<= 41.7.
train96catBSEP<-train96p$BSEP.Inhibitor
names(train96catBSEP)<-train96p$SampleID

##To map between expression and phenotype files.
##colnames(EXPRESSION VARIABLES) = rownames(PHENOTYPE AND METADATA FILES) = names(train96IC50, train96catBSEP)

  
#########################################
######TRAINING DATA FOR 24-HOUR EXPOSURES
#########################################
#Contains 237 drug response signatures, representing 79 compounds (9 BSEP inhibitors and 65 noninhibitors)

#Expression response to drug exposure
train24<-loadEntity(qryTRAIN[6,2])
train24fc<-train24$objects$FoldChangeResponse  #Drug response signature calculated as fold change (H/aveC)
train24ld<-train24$objects$LogDeltaResponseData #Drug response signature calculated as difference in logs (logH - logaveC)
probeMap<-train24$objects$ProbeAnnotations
#Phenotypes
trainp24<-loadEntity(qryTRAIN[4,2])
train24p<-trainp24$objects$Phenotypes  #BSEP phenotypes
train24m<-trainp24$objects$Metadata  #information about originating CEL files and the rats.

#Continuous trait for BSEP inhibition (IC50 for rat BSEP inhibition using in vitro assay)
train24IC50<-train24p$RatBSEPIC50
names(train24IC50)<-train24p$SampleID
#Categorical trait for BSEP inhibition.  0 where ratIC50>=100 and 1 where ratIC50<= 41.7.
train24catBSEP<-train24p$BSEP.Inhibitor
names(train24catBSEP)<-train24p$SampleID

##To map between expression and phenotype files.
##colnames(EXPRESSION VARIABLES) = rownames(PHENOTYPE AND METADATA FILES) = names(train24IC50, train24catBSEP)


########################TRANING DATA####################
##Test data is located in folder: syn1468168
##Data was provided by Amgen Discovery Toxicology Department
##Data was normalized using SNM and Drug Response signatures were calculated separately for compound 
#by calculating either the fold change (H/aveC) or Logged Difference (logH - logaveC)
#for each expression signature following high dose exposure (H) relative the average expression for the three control-exposed samples (aveC)

qryTEST <- synapseQuery('select id, name from entity where entity.parentId=="syn1468168"')

#########################################
######TEST DATA FOR 96-HOUR and 72-HOUR EXPOSURES
#########################################
#Contains 75 drug response signatures, representing 17 compounds (7 BSEP inhibitors and 9 noninhibitors)
#For each compound, data was collected following 72-hour (for twelve of the compounds) or 96-hours (for five of the compounds)  in multiple rats each for exposure to a 'high' dose of compound or control.  
#Expression profiles were measured using the Affy rat 230_2 chip (which is the same used in the training data) expect for 3 of the 72-hour exposures, which used the Affy rat 230a chip that contains a subset of ~80% of the probes represented on the 230_2 chip.

#Expression response to drug exposure
test96<-loadEntity(qryTEST[1,2])
test96fc<-test96$objects$FoldChangeResponse  #Drug response signature calculated as fold change (H/aveC)
test96ld<-test96$objects$LogDeltaResponseData #Drug response signature calculated as difference in logs (logH - logaveC)
probeMap<-test96$objects$ProbeAnnotations
#Phenotypes
testp96<-loadEntity(qryTEST[5,2])
test96p<-testp96$objects$Phenotypes  #BSEP phenotypes
test96m<-testp96$objects$Metadata  #information about originating CEL files and the rats.

#Continuous trait for BSEP inhibition (IC50 for rat BSEP inhibition using in vitro assay)
test96IC50<-test96p$RatBSEPIC50
names(test96IC50)<-test96p$SampleID
#Categorical trait for BSEP inhibition.  0 where ratIC50>=100 and 1 where ratIC50<= 41.7.
test96catBSEP<-test96p$BSEP.Inhibitor
names(test96catBSEP)<-test96p$SampleID

##To map between expression and phenotype files.
##colnames(EXPRESSION VARIABLES) = rownames(PHENOTYPE AND METADATA FILES) = names(test96IC50, test96catBSEP)

#########################################
######TEST DATA FOR 24-HOUR EXPOSURES
#########################################
#Contains 142 drug response signatures, representing 37 compounds (10 BSEP inhibitors and 22 noninhibitors)
#For each compound, data was collected following 24-hours in 2-5 rats each for exposure to a 'high' dose of compound or control.  Expression profiles were measured using the Affy rat 230_2 chip (which is the same used in the training data) for all expect for 3 of the compounds, which used the Affy rat 230a chip that contains a subset of ~80% of the probes represented on the 230_2 chip.

#Expression response to drug exposure
test24<-loadEntity(qryTEST[6,2])
test24fc<-test24$objects$FoldChangeResponse  #Drug response signature calculated as fold change (H/aveC)
test24ld<-test24$objects$LogDeltaResponseData #Drug response signature calculated as difference in logs (logH - logaveC)
probeMap<-test24$objects$ProbeAnnotations
#Phenotypes
testp24<-loadEntity(qryTEST[4,2])
test24p<-testp24$objects$Phenotypes  #BSEP phenotypes
test24m<-testp24$objects$Metadata  #information about originating CEL files and the rats.

#Continuous trait for BSEP inhibition (IC50 for rat BSEP inhibition using in vitro assay)
test24IC50<-test24p$RatBSEPIC50
names(test24IC50)<-test24p$SampleID
#Categorical trait for BSEP inhibition.  0 where ratIC50>=100 and 1 where ratIC50<= 41.7.
test24catBSEP<-test24p$BSEP.Inhibitor
names(test24catBSEP)<-test24p$SampleID

##To map between expression and phenotype files.
##colnames(EXPRESSION VARIABLES) = rownames(PHENOTYPE AND METADATA FILES) = names(test24IC50, test24catBSEP)
