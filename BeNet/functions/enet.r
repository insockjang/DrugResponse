lambda.enet.cv <- function(X, y, K, l1range, l2range){
	
	lambda.list <- matrix(0, nrow=length(l2range), ncol=3)
	
	for(i in 1:length(l2range)){
		enet.cv <-  cv.enet(X, y, K, plot.it = FALSE, mode = "fraction", s = l1range, lambda = l2range[i], intercept=FALSE, normalize=FALSE);
		lambda.list[i,] <- c(enet.cv$s[enet.cv$cv == min(enet.cv$cv)], l2range[i], min(enet.cv$cv))
	}
	
	lambda <- lambda.list[lambda.list[,3] == min(lambda.list[,3]),]
	
	lambda[1:2] # return lambda[1] and lambda[2]	
}


calc.reg.coef.enet <- function(X, y, lambda){
	
	model <- enet(X, y, intercept = FALSE, lambda[2], normalize = FALSE)
	coef <- predict.enet(model, type = "coefficients", mode = "fraction", s = lambda[1])$coefficients
	
	as.matrix(coef) # return the coefficients beta
}