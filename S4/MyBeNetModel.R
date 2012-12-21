# Warning : Do not run this.
# I will post local function repo (source code) for implementing Bayesian elastic Net customized class 
# Need to modify
# The REAL customized(wrapped up) class $ fuction outside caret
# Bayesian Elastic Net Regression Model is included into our predictiveModeling
# enetqr function generate class which is 'list'

setClass(
  Class = 'MyBeNetModel',
  contains="PredictiveModel",
  representation = representation(beNetModel='list'))


# we only need beta and alpha(intercepts) from MCMC samples
# We implemented 1000 MCMC samples and discarded Burn-in period till first 200 samples.
# This is for custom training

setMethod(

                f = "customTrain",

                signature = signature("MyBeNetModel", "matrix", "numeric"),

                definition = function(method, featureData, responseData){

                        message("in train_BayesianElasticNetRegression")
                        
                        beNetModel <- model.selection(responseData ~ featureData,tau=0.5,burnin = 2, mcmc = 10)
                        
                        coeff <- beNetModel$Top.model
                        beta <- coeff[-1]
                        names(beta) <- colnames(featureData)
          
                        beNetModel$coef <- beta
                        beNetModel$intercept <- coeff[1]
                      
                        method@beNetModel <- beNetModel

                        method

                }

)

# This is for custom prediction
# Used with estimated coefficients(including beta and intercept) from custom training

setMethod(

                f = "customPredict",

                signature = signature("MyBeNetModel", "matrix"),

                definition = function(method, featureData){

                        print("in predict_BayesianElasticNetRegresion")

                        testCoefs <- method@beNetModel$coef

                        testCoefs <- as.matrix(testCoefs)

                        testFeatures <- featureData[, row.names(testCoefs)]

                        prediction <- testFeatures %*% testCoefs

                        prediction <- prediction + method@beNetModel$intercept

                        return(prediction)

                }

)    

 

MyBeNetModel <- function(){

 newObj <- new("MyBeNetModel")

}

 