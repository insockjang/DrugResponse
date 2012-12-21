# Sample code Curation with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
## load the synapse client and login
library(synapseClient)
library(affy)
synapseLogin()

## set up a project
myName <- "Sanger_Exp_GCRMA"
projName <- sprintf("%ss Curation Project %s", myName, as.character(gsub("-",".",Sys.Date())))
myProj <- createEntity(Project(list(name=projName)))

# This is parent ID(project id)
Sanger_parentID <- "5019"

## create a dataset 
myDataset <- createEntity(Dataset(list(name="Sanger Expression GCRMA", parentId=Sanger_parentID)))

## create a layer if parent id is mapped into a dataset
#myLayer <- createEntity(Layer(list(name="Sanger Expression RMA", type="E", parentId=Sanger_parentID)))
## view the dataset on the web to add a description
onWeb(myDataset)

## refresh the local copy of myDataset
myDataset <- refreshEntity(myDataset)


################
## Section 2: Working with data
################

## create a new expression layer
myExpr <- createEntity(Layer(list(name="curated expression", type="E", parentId = propertyValue(myDataset, "id"), status="curated")))

## add an annotation specifying the data format
annotValue(myExpr, "format") <- "sageBioCurated_withROBJECT"

## add the pm data file to the entity
myExpr <- addObject(myExpr, FIT.RMA)

## store the data
myExpr <- storeEntity(myExpr)

## RObject upload
probeIDLayer <- Layer(list(name = "R_Expression_Sanger_GCRMA", type = "E", parentId = "104242", status="db"))
testLayer <- addObject(probeIDLayer, FIT.GCRMA)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Expression_Sanger_RMA", type = "E", parentId = "104230", status="db"))
testLayer <- addObject(probeIDLayer, FIT.RMA)
testLayer <- storeEntity(testLayer)

probeIDLayer <- Layer(list(name = "R_Expression_Sanger_metaGEO", type = "E", parentId = "104395", status="db"))
testLayer <- addObject(probeIDLayer, fit)
testLayer <- storeEntity(testLayer)

