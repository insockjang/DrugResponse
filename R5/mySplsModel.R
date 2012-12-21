#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
require(spls)
mySplsModel <- setRefClass(Class = "mySplsModel",
                           contains="PredictiveModel",
                           fields="model",
                           methods = list(
                             initialize = function(...){
                               return(.self)
                             },
                             
                             rawModel = function(){
                               return(.self$model)
                             },
                             
                             customTrain = function(featureData, responseData, nfolds = 5, K, eta, ...){
                               MAT<-matrix(0,nrow=length(eta),ncol=length(K))
                               for(k1 in 1:length(eta)){
                                 for(k2 in 1:length(K)){
                                   foldIndices <- createFolds(featureData[,1], k = nfolds, list = TRUE)                                   
                                   
                                   foldPredict <- foreach(fold = foldIndices) %dopar% {
                                     filterData<-filterPredictiveModelData(featureData[-fold,],responseData[-fold])
                                     featureData1<-filterData$featureData
                                     responseData1<-filterData$responseData
                                     fit<-spls(featureData1,responseData1,eta=eta[k1],K=K[k2])
                                     testPredict<-predict(fit,data=featureData[fold,])
                                     MSE = mean((testPredict - responseData[fold])^2)
                                     return(MSE)
                                   }
                                   
                                   MAT[k1,k2]<-mean(do.call("c",foldPredict))
                                 }
                               }
                               rownames(MAT)<-eta
                               colnames(MAT)<-K
                               pos<-which(MAT == min(MAT), arr.ind = TRUE)
                               eta_opt<-eta[pos[1]]
                               K_opt<-K[pos[2]]
                               
                               filterData1<-filterPredictiveModelData(featureData,responseData)
                               featureData2<-filterData1$featureData
                               responseData2<-filterData1$responseData
                               
                               .self$model<-spls(featureData2,responseData2,eta=eta_opt,K=K_opt)                                 
                             },
                             
                             customPredict = function(featureData){
                               predictedResponse <- predict(.self$model, featureData)
                               return(predictedResponse)
                             },
                             
                             getCoefficients = function(){
                               coeff<-predict(.self$model,type = "coefficient")
                               return(coeff)
                             }
                             
                             )
                           
                           )
