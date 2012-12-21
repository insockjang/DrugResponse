# for testing customized Class for least angle regression as one testbed for predictiveModeling

setClass(
  Class = 'MyLarsModel',
  contains="PredictiveModel",
  representation = representation(larsModel='lars')
)

setMethod(

                f = "customTrain",

                signature = signature("MyLarsModel", "matrix", "numeric"),

                definition = function(method, featureData, responseData){

                        message("in customTrain ", class(method))

                        larsModel <- lars(featureData,responseData,type = "lar",use.Gram = FALSE, max.steps=floor(ncol(featureData) * 0.1)
                        
                        beta <- larsModel$beta[nrow(larsModel$beta),]
                        
                        names(beta) <- colnames(featureData)
          
                        intercept <- larsModel$mu
                        print(intercept)
                        print(length(beta))

                        larsModel$coef <- beta
                        larsModel$intercept <- intercept
                      
                        method@larsModel <- larsModel

                        method

                }

)

# This is for custom prediction
# Used with estimated coefficients(including beta and intercept) from custom training

setMethod(

                f = "customPredict",

                signature = signature("MyLarsModel", "matrix"),

                definition = function(method, featureData){

                        message("in customPredict ", class(method))

                        testCoefs <- method@larsModel$coef

                        testCoefs <- as.matrix(testCoefs)

                        testFeatures <- featureData[, row.names(testCoefs)]

                        print(dim(testFeatures))
                        
                        prediction <- testFeatures %*% testCoefs

                        prediction <- prediction + method@larsModel$intercept

                        return(prediction)

                }

)    

 

MyLarsModel <- function(){

 newObj <- new("MyLarsModel")

}
