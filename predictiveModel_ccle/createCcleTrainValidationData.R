### DEMO SVM
library(predictiveModeling)
library(synapseClient)
source("/home/ijang/COMPBIO/trunk/users/jang/predictiveModel_ccle/createCcleTrainValidationData.R")
###################################################
#### Load Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "48339"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$copySet

id_oncomapLayer <- "48341"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$oncomapSet

id_exprLayer <- "48344" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$exprSet

###################################################
### Load Response Data
###################################################

id_drugLayer <- "48359" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$responseADF

###### how many times to validate implemented model 
numVals = 5


######### store all competition data in Synapse ####################
project_competition <- Project(list(name = "CCLE predictiveModeling Comparison data", description = "This project contains all data with 5 randomly chosen"))
project_competition <- createEntity(project_competition)

dataset_ccle <- Dataset(list(name = "CCLE dataset with train and validate", parentId = propertyValue(project_competition, "id")))
dataset_ccle <- createEntity(dataset_ccle)

for(k in 1:numVals){
  trValDataSet <- createCcleTrainValidationData(exprData = eSet_expr, copyData = eSet_copy, oncomapData = eSet_oncomap, drugData = adf_drug)  
  
  data <- Layer(list(name = paste("CCLE dataset_0",k,sep = ""), type = "E", parentId = propertyValue(dataset_ccle, "id")))
  data <- addObject(data, trValDataSet)
  data <- storeEntity(data)
  
}
