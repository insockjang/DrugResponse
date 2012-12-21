# Still to do: ####### Initialization and Libraries ####### 
####### ---------------------------- ####### 
rm(list=ls()); 
library(lars);		#LASSO 
library(elasticnet);	#Elastic Net, requires package "lars" 
library(MASS); 
library(statmod); 
library(MCMCpack); 
library(predictiveModeling)
#Read in functions for 
for (nm in list.files("./functions", pattern = "\\.[RrSsQq]$")) {          
	source(file.path("./functions", nm))
}
	
####### Import Data ############ 
####### ----------- ############ 

# Read in data: y is an "n x 1" vector of responses centered at 0, 
# X is an "n x p" design matrix with no intercept column 


library(lasso2) 

data(demoData)


ds_features_cn_mut_ccle <- createAggregateFeatureDataSet(list(copy = copySet, mut = oncomapSet))
checkEquals(nrow(ds_features_cn_mut_ccle), nrow(exprs(copySet)) + nrow(exprs(oncomapSet)))

dataSets_ccleFeaturesCnMut_sangerChems <- createFeatureAndResponseDataList(ds_features_cn_mut_ccle, sangerADF)

DATA <- filterPredictiveModelData(t(dataSets_ccleFeaturesCnMut_sangerChems$featureData),dataSets_ccleFeaturesCnMut_sangerChems$responseData["PLX4720",,drop=FALSE])



yfull = as.matrix(DATA$responseData); 
Xfull = as.matrix(DATA$featureData); 

num.of.sets = 10;
MSE = mat.or.vec(num.of.sets,7)
Betas = array(0,dim=c(num.of.sets,7,dim(Xfull)[2]))

NUM<-floor(length(yfull)*0.7)

for(l in 1:num.of.sets){  #100 different splits of the data

index = sample(1:length(yfull))

X <- Xfull[index[1:NUM],] 
y <- yfull[index[1:NUM]] 

y<-scale(y)
X <-scale(X)

Xtest <- Xfull[index[NUM+1:length(yfull)],]; 
ytest <- yfull[index[NUM+1:length(yfull)]]; 

ytest <- scale(ytest)
Xtest <- scale(Xtest)
	

####### Elastic Net (classical) ##########
####### ----------------------- ##########

# Range to search over using CV
l1range = (0:100)/100;
l2range = 2^((-10):10);

#Calculate lambdas with lowest CV error
lambda.enet = lambda.enet.cv(X, y, K=5, l1range, l2range);

#Calculate regression coefficients
Betas.Enet = calc.reg.coef.enet(X, y, lambda.enet);

####### Bayesian Elastic Net ##########
####### -------------------- ##########

lambda1.bayesenet = max(abs(t(X) %*% (y-X %*% Betas.Enet)));				#Use LARS MAP estimate to find lambda1, lambda2 with lowest CV error
lambda2.bayesenet = lambda.enet[2] + .001;	

#Calculate regression coefficients
Betas.Bayesenet = calc.reg.coef.bayesenet(X, y, lambda1.bayesenet, lambda2.bayesenet, 10000);

######## Plots ############
######## ----- ############

#plot(Betas.LS)
#points(Betas.Normal, col="blue")
#points(Betas.Lasso, col="red")
#points(Betas.Enet, col="green")
#points(Betas.Bayeslasso, col="purple")
#points(Betas.Bayesenet, col="yellow")
#points(Betas.Bayeslq, col="lightblue")

Betas[l,1,] = Betas.Enet;
Betas[l,2,] = Betas.Bayesenet;

MSE[l,1] = mean((ytest - Xtest%*%Betas.Enet)^2)
MSE[l,2] = mean((ytest - Xtest%*%Betas.Bayesenet)^2)

save(MSE, Betas, file="MSE_Betas.RData")

}