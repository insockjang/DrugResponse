#' Constructor for a class in the PredictiveModel class hierarchy that implements a simple function regressing the top correlated features
#' R5 Class for Bayesian Lasso for predictiveModeling
#' Bayesian Lasso Regression Model is now included into our predictiveModeling
#'
#' @export
require("standGL")
require("synapseClinet")
myOverlapGroupLassoModel <- setRefClass(Class = "myOverlapGroupLassoModel",
                                        contains="PredictiveModel",
                                        fields=list(coefficients="numeric"),
                                        methods = list(
                                          initialize = function(...){
                                            return(.self)
                                          },
                                          
                                          train = function(featureData, responseData, trControl = defaultTrainControl(), filterData = TRUE, tuneGrid = NULL,index=NULL, ...){
                                            if(filterData == TRUE){
                                              message("filtering data...")
                                              processedData <- filterPredictiveModelData(featureData, responseData)
                                              featureData <- processedData$featureData
                                              responseData <- processedData$responseData
                                            }
                                            
                                            grouping = function(featureData){
                                              entity<-loadEntity(105350)
                                              mSigDB_annotations <- loadEntity(105363)                                              
                                              DB<-entity$objects$MsigDB_symbolID                                              
                                              allPathways <- mSigDB_annotations$objects$C2$KEGG
                                              DATA<-featureData
                                              DATA1<-c()
                                                                                        
                                              index <-c()
                                              j<-0
                                              for (i in 1:length(allPathways)){
                                                curPathway1 <- allPathways[[i]]
                                                mSigDB_index1 <- which(DB$geneset.names == curPathway1)
                                                curPathwayGenes1 <- DB$genesets[[mSigDB_index1]]
                                                AA<-rownames(DATA) %in% curPathwayGenes1
                                                pathwayGene <- DATA[which(AA==1),]
                                                if (length(dim(pathwayGene)[1]) !=0){
                                                  j<-j+1
                                                  a<-dim(pathwayGene)[1]  
                                                  b<-rep(j,a)
                                                  index <-c(index,b)
                                                  DATA1 <-rbind(DATA1,pathwayGene)
                                                }
                                              }                                             
                                              returen(list(groupData=DATA1,index = index)))
                                            }
                                            
                                            # add intercept part
                                            DATA <- rbind(rep(1,dim(grouping$groupData)[2]),grouping$groupData)
                                            
                                            # kind of tuneGrid for lambda
                                            
                                            fit <- standGL(responseData,DATA,index=index,alpha = 1,is.pen=as.numeric(c(0,index>0)),lam.path = tuneGrid)
                                            
                                            Y <-fit$beta %*% DATA
                                            
                                            MSE<-c()
                                            for (i in 1:dim(Y)[2]){
                                              MSE<-c(MSE,mean((Y[,i]-responseData)^2))
                                            }
                                            
                                            .self$lambdaOpt <- fit$lam.path[which(MSE == min(MSE))]                                 
                                            .self$coefficients <- fit$beta[,which(MSE == min(MSE))]
                                          },
                                          
                                          
                                          predict = function(featureData){
                                            # Be careful when p>>n : This is singular case. So we have to consider "NA" for predicting response
                                            featureNames <- names(.self$coefficients[2:length(.self$coefficients)])
                                            coefNames <-names(.self$coefficients[which(.self$coefficients != 0)])
                                            interNames <- intersect(featureNames,coefNames)
                                            predictedResponse <- .self$coefficients[1] + featureData[,interNames] %*% .self$coefficients[interNames]
                                            return(predictedResponse)
                                          }
                                          )
                                        )
