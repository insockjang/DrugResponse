Bridge <- function(B,y,X,lam,gam) {   
  ## Rosenbrock Banana function
 t(y - X %*% B) %*% (y - X %*% B) + lam * sum(abs(B)^gam)
}


bridge.CV <- function(X, y, lam.range, gam.range, folds = dim(X)[1]){

    n <- dim(X)[1] 
    p <- dim(X)[2]
    error <- array(0, dim<-c(folds,length(lam.range),length(gam.range)))

	indices <- sample(1:n,folds)
    
    for(i in 1:length(lam.range)){
        lam <- lam.range[i];
        for(j in 1:length(gam.range)){
            gam <- gam.range[j]
            for(k in 1:folds){
					      k.ind <- indices[k]
                Xtrain <- X[-k.ind,] 
                ytrain <- y[-k.ind]
                Xtest <- X[k.ind,] 
                ytest <- y[k.ind]
                
                reg.coef <- optim(Betas.LS, Bridge, NULL, ytrain, Xtrain, lam, gam)$par;
                error[k,i,j] <- (t(reg.coef) %*% Xtest - ytest)^2;
            }
        }
    }
    error.mean <- apply(error,2:3, mean);
    error.min1 <- apply(error.mean, 1, min);
    error.min2 <- apply(error.mean, 2, min);
    lam.final <- lam.range[which(error.min1 %in% min(error.min1))];
    gam.final <- gam.range[which(error.min2 %in% min(error.min2))];  

    c(lam.final, gam.final) ## return
}


## Gussian Prior
gPrior.cv <- function(X, y, g.range){

    n <- dim(X)[1]
    p <-dim(X)[2]
    error <- mat.or.vec(n,length(g.range))
    
    for(i in 1:length(g.range)){
        g <- g.range[i]
        for(j in 1:n){
            Xtrain <- X[-j,]
            ytrain <- y[-j]
            
            Xtest <- X[j,]
            ytest <- y[j]
            
            Sigma <- g * solve(t(Xtrain) %*% Xtrain)
            reg.coef <- solve(Sigma + t(Xtrain) %*% Xtrain) %*% (t(Xtrain) %*% ytrain)
            error[j,i] <- (t(reg.coef) %*% Xtest - ytest)^2
        }
	}
    error.mean <- apply(error,2,mean)
    g <- g.range[which(error.mean%in%min(error.mean))]
}


calc.reg.coef.bayeslq <- function(X, y, lambda, q, M){

  ##### Constants and variables #####
  ##### ----------------------- #####
  n <- length(y);
  p <- ncol(X);	
  M <- M;		#length of Markov chain

  a <- 1
  b <- 1  #Inverse Gamma prior parameters for sigma^2


  betas <- matrix(0, nrow<-M,ncol<-p);
  sigma2 <- rep(0,M);
  lam1 <- lambda;
  q <- q;
  	
  ##### MCMC Time 1 ######
  ##### ------ ######
  
  sigma2[1] <- var(y)
  betas[1,] <- solve(t(X)%*%X)%*%t(X)%*%y
  
  ##### MCMC Time<-2,...,M #####
  ##### ----------------- #####
  
  
  #Update Parameters
  for (i in 2:M) {
    #update betas
    mean <- solve(t(X) %*% X) %*% t(X) %*% y
    sig1 <- sigma2[i-1] * solve(t(X) %*% X)
    betasprop <- mvrnorm(1, mean, sig1)
          
    acc <- exp((-lam1)/(2*sigma2[i-1]) * (sum(abs(betasprop)^q) - sum(abs(betas[i-1,])^q) ) )
          
    if(runif(1)< acc){
      betas[i,] <- betasprop
    }
    else{
      betas[i,] <- betas[i-1,]
    }
          
          
  
    #update sigma
    sigshape <- (n+p)/2 + a
    sigscale <- t(y-X%*%betas[i,])%*%(y-X%*%betas[i,])/2 + (lam1/2)*sum(abs(betas[i,])^q)  + b
    sigma2[i] <- rinvgamma(1, shape <- sigshape, scale <- sigscale) 
  
  }
  
  	
  coef <- apply(betas, 2, mean)
  coef
    
}
