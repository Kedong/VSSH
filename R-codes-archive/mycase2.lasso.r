rm(list=ls())

setwd('/Users/swang/Box Sync/Academic/Research/Projects/William Li/heredity/interaction')
library(glmnet)

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

q = 10
beta1 = c(1,1,1,1,rep(0,6))/10
beta2 = rep(0,45)
beta2[c(1,2,3,10,11,18)] = rep(10,6)
beta = c(beta1, beta2)
sigma = 12

for (sim in 1:100) {

	cat('sim:', sim, '\n')

	set.seed(sim*100+12345)

	n = 200

	## generate training data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)
	x.tr = x
	y.tr = y

	## generate validation data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)
	x.val = x
	y.val = y

	## generate test data
	n = 10000
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)
	x.te = x
	y.te = y

	fit = glmnet(x=x.tr, y=y.tr)
	yhat.val = predict(fit, x.val)
	idx = order(apply((yhat.val-as.vector(y.val))^2,2,mean))[1]

	ahat = as.numeric(fit$a0[idx])
	bhat = as.vector(fit$beta[,idx])

	yhat.te = ahat + x.te%*%bhat
	mse.te = mean((yhat.te-y.te)^2)

	write.table(matrix(c(ahat,bhat),nrow=1),'coefmat.lasso.mycase2.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
	write.table(mse.te,'infomat.lasso.mycase2.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
}
