gPrior.cv <- function(X, y, g.range){

  n <- dim(X)[1]
  p <- dim(X)[2];
  error <- mat.or.vec(n,length(g.range))
    
  for(i in 1:length(g.range)){
    g <- g.range[i]
    for(j in 1:n){    
      Xtrain <- X[-j,]
      ytrain <- y[-j]
      Xtest <- X[j,]
      ytest <- y[j]
            
      Sigma <- g * solve(t(Xtrain) %*% Xtrain)
      reg.coef <- solve(Sigma + t(Xtrain)%*%Xtrain) %*% (t(Xtrain) %*% ytrain)
      error[j,i] <- (t(reg.coef) %*% Xtest - ytest)^2
    }
  }
  error.mean <- apply(error,2,mean)
  g <- g.range[which(error.mean %in% min(error.mean))]
}


calc.reg.coef.normal <- function(X, y, Sigma){
	
	#reg.coef = solve(Sigma + t(X)%*%X) %*% ((t(X)%*%X)%*%(solve(t(X)%*%X)%*%t(X)%*%y))

	reg.coef <- solve(Sigma + t(X)%*%X) %*% (t(X)%*%y)
	
	reg.coef # return beta coefficients
}
