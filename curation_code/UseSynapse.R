# Create Project 
library(synapseClient)
synapseLogin()

## set up a project
myName <- "Sanger Data"
projName <- sprintf("%ss Curation Sanger %s", myName, as.character(gsub("-",".",Sys.Date())))

myProj <- createEntity(Project(list(name=projName)))

# This is parent ID(project id)
Sanger_parentID <- "5019"

## create a dataset 
myDataset <- createEntity(Dataset(list(name="Sanger Dataset", parentId=Sanger_parentID)))

## create a layer if parent id is mapped into a dataset
#myLayer <- createEntity(Layer(list(name="Sanger Expression RMA", type="E", parentId=Sanger_parentID)))

## view the dataset on the web to add a description
onWeb(myDataset)

## refresh the local copy of myDataset
myDataset <- refreshEntity(myDataset)

#################### GSK
myName <- "GSK Data"
projName <- sprintf("%ss Curation GSK %s", myName, as.character(gsub("-",".",Sys.Date())))

myProj <- createEntity(Project(list(name=projName)))

# This is parent ID(project id)
GSK_parentID <- "5019"

## create a dataset 
myDataset <- createEntity(Dataset(list(name="GSK Dataset", parentId=GSK_parentID)))

## create a layer if parent id is mapped into a dataset
#myLayer <- createEntity(Layer(list(name="Sanger Expression RMA", type="E", parentId=Sanger_parentID)))

## view the dataset on the web to add a description
onWeb(myDataset)

## refresh the local copy of myDataset
myDataset <- refreshEntity(myDataset)




# syn377397
load("~/PredictiveModel/drugResponse/CCLE/NP/bootstrapENet2_woUnique/reference.Rdata")

load("~/PredictiveModel/drugResponse/CCLE/NP/bootstrapENet2_hybrid_woUnique/reference.Rdata")

projName <- "Cancer Cell Line Project - Molecular Features"

myProj <- createEntity(Project(list(name=projName)))

myFolder <- createEntity(Folder(list(name = "CCLE", parentId = "syn1528791")))

myFolder1<-createEntity(Folder(list(name = "EXP + COPY + MUT(oncomap)",parentId = myFolder$properties$id)))

newEntity <- Data(list(name= "Selected Molecular Features", parentId = myFolder1$properties$id))

newEntity<-createEntity(newEntity)          

addObject(newEntity, reference)

storeEntity(newEntity)

myFolder1<-createEntity(Folder(list(name = "EXP + COPY + MUT(hybrid)",parentId = myFolder$properties$id)))

newEntity <- Data(list(name= "Selected Molecular Features", parentId = myFolder1$properties$id))

newEntity<-createEntity(newEntity)          

addObject(newEntity, reference)

storeEntity(newEntity)