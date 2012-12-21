calc.reg.coef.bayesenet <- function(X, y, lambda1, lambda2, M){
  
  ##### Constants and variables #####
  ##### ----------------------- #####
  n <- length(y)
  p <- ncol(X)
  M <- M		#length of Markov chain
  
  a <- 1
  b <- 1  #Inverse Gamma prior parameters for sigma^2
  
  
  betas <- matrix(0, nrow = M, ncol = p)
  taus <- matrix(0, nrow = M, ncol = p)
  sigma2 <- rep(0,M)
  lam1 <- lambda1
  lam2 <- lambda2
  	
  ##### Time 1 ######
  ##### ------ ######
  
  sigma2[1] <- var(y)
  taus[1,] <- rep(1,p)
  betas[1,] <- solve(t(X) %*% X) %*% t(X) %*% y
  
  ##### MCMC Time=2,...,M #####
  ##### ----------------- #####
  
  
  #Update Parameters
  for (i in 2:M) {
    D <- diag(1/taus[i-1,],ncol=p)
    A <- t(X)%*%X + D + diag(rep(lam2,p))
    Ainv <- solve(A)

    #update betas
    mean <- Ainv %*% t(X) %*% y
    sig1 <- sigma2[i-1] * Ainv
    betas[i,] <- mvrnorm(1, mean, sig1)

    #update sigma
    sigshape <- (n-1+p)/2 + a
    sigscale <- t(y-X%*%betas[i,])%*%(y-X%*%betas[i,])/2 + t(betas[i,])%*%D%*%betas[i,]/2 + t(betas[i,])%*%diag(rep(lam2,p))%*%betas[i,]/2 + b
    sigma2[i] <- rinvgamma(1, shape = sigshape, scale = sigscale) 

    #update taus
    for(j in 1:p){
      taus[i,j] <- 1 / rinvgauss(1, mu = sqrt(lam1^2*sigma2[i]/(betas[i,j]^2)), lambda = lam1^2)
    }
  }
	
	coef <- apply(betas, 2, mean);
	coef #return beta

}