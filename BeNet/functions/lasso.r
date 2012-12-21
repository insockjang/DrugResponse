# library lars contains "lasso"
# using cross validation to compute lambda
lambda.lasso.cv <- function(X, y, K){
	
	lasso.cv <-  cv.lars(X, y, K, plot.it = FALSE, intercept = FALSE, normalize = FALSE);
	
	lambda <- lasso.cv$fraction[lasso.cv$cv == min(lasso.cv$cv)];
	
	lambda;	
}

# given lambda, find beta coefficients 
calc.reg.coef.lasso <- function(X, y, lambda){
	
	model <- lars(X, y, type<- "lasso", intercept = FALSE, normalize = FALSE);
	coef <- coef.lars(model, mode = "fraction", s = lambda);
	
	as.matrix(coef);
}
