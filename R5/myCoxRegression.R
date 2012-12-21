#' Constructor for a class in the PredictiveModel class hierarchy that wraps models returned from caret
#'
#' @param model a raw model returned by caret train()
#' @param modelType
#' @return an instance of the class CaretModel
#' @seealso caret
#' @export
require(glmnet)

myPenaltyCoxModel <- setRefClass(Class = "myPenaltyCoxModel",
                        contains="PredictiveModel",
                        fields="model",
                        methods = list(
                             initialize = function(...){
                               return(.self)
                             },
                             
                             rawCaretModel = function(){
                               return(.self$model)
                             },
                             
                             train = function(featureData, responseData, trControl = defaultTrainControl(),
                                              filterData = TRUE, tuneGrid = NULL, penaltyFactor = rep(1,dim(featureData)[2])){
                                .self$model <- cv.glmnet(featureData,responseData,family = "cox",alpha = 1,nfolds = 3, penalty.factor = penaltyFactor) 
                                
                             },
                             
                             customPredict = function(featureData){
                               predictedResponse <- predict(fit,featureData,s="lambda.min")
                               return(predictedResponse)
                             }
                             
                             
                             )
                           )
