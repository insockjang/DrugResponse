bootstrapFeatureSelection <- 
  function(model, featureData, responseData, numBootstrap = 100,  ...){
        
    bootIndices <- createResample(featureData[,1],times = numBootstrap,list = TRUE)
    
    bootResults <- foreach (numBoot = bootIndices) %dopar% {
      bootModel<-model$new()
      bootModel$customTrain(featureData[numBoot, ], responseData[numBoot], ...)
      coeff <- bootModel$getCoefficients()      
      return(as.matrix(coeff!=0))
    }
    
    sumBootCoef<-matrix(0, nrow = ncol(featureData),ncol = ncol(bootResults[[1]]))
    for(k in 1:length(bootResults)){
      sumBootCoef<-sumBootCoef+bootResults[[k]]
    }
    return(sumBootCoef)
  }

bootstrapFeatureGrid <-
  function(model, featureData, responseData, numBootstrap = 100, alpha=1, ...){
    
    alphaBootstrap <- foreach(k = 1:length(alpha)) %dopar% {
      bootCoeff<-bootstrapFeatureSelection(model, featureData, responseData, numBootstrap = 100, alpha=alpha[k],  ...)
      return(bootCoeff)      
    }
    return(alphaBootstrap)
  }
