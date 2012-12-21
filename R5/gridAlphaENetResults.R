gridAlphaENetResults<-function(featureData, responseData, model, alphas = alphas, lambda = lambda, ...){
  
  alphaResults <- foreach(k1 = 1:length(alphas)) %dopar% {
    alpha <- alphas[k1]
    foldModel <- model$new()
    foldModel$customTrain(featureData, responseData, alpha = alpha, lambda = lambda, ...)
    return(foldModel)
  }
  
  coeffs<-c()
  cvms <-c()
  alpha_temp <-c()
  lambda_temp <-c()
  for(i in 1:length(alphas)){
    fit<-alphaResults[[i]]$rawModel()
    MEDIAN <- fit$cvm
    names(MEDIAN) <- fit$lambda
    
    alpha_temp <- c(alpha_temp,fit$cvm[which.min(fit$cvm)])
    lambda_temp <- c(lambda_temp,fit$lambda[which.min(fit$cvm)])
        
    cvms <- rbind(cvms,MEDIAN)    
    coeffs <- cbind(coeffs,alphaResults[[i]]$getCoefficients())    
  }   
  
  optAlpha <- alphas[which.min(alpha_temp)]
  optLambda <- lambda_temp[which.min(alpha_temp)]
  
  rownames(cvms) <- alphas
  colnames(coeffs) <- alphas
  return(list(cvm = cvms, optAlpha = optAlpha, optLambda = optLambda, bestModel = alphaResults[[which.min(alpha_temp)]]))
    
}
