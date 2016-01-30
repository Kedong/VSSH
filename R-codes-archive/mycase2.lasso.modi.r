rm(list=ls())

setwd('D:/Schooling/Higher Education/!papers&research/SHIM_New')
library(glmnet)

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numintzero = 0  # number of interactions that were set as 0
numnonusemain = 0  # number of mains that dont have interactions
constcoef = 0  # if the coefficients are constant
main.magnitude = 10
int.magnitude = 10

beta1 = c(rep(1,numsig),rep(0,q-numsig))
beta2 = combn(c(rep(1,numsig-numnonusemain),rep(0,q-numsig+numnonusemain)),2,prod)
if (numintzero!=0) beta2[sample(which(beta2>0),numintzero)]=0
beta1 = beta1 * main.magnitude
beta2 = beta2 * int.magnitude
beta = c(beta1, beta2)
if (constcoef!=1) {for (i in 1:length(beta)) beta[i]=beta[i]*(runif(1)+0.5)}
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

	write.table(matrix(c(ahat,bhat),nrow=1),'coeflasso.randcoef.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
	write.table(mse.te,'infolasso.randcoef.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
}
