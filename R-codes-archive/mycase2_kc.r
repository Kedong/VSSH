rm(list=ls())

setwd('D:/Schooling/Higher Education/!papers&research/SHIM_New')
library(glmnet)
source('selection_interaction_kc.r')

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

q = 10  # 10 factors, so 10C2 = 45 2-factor-interactions
beta1 = c(1,1,1,1,rep(0,6))/10 # main.coef = 0.1
beta2 = rep(0,45)
beta2[c(1,2,3,10,11,18)] = rep(10,6) # the CORRESPONDING inter.coef = 10, see if LASSO can pick them out
beta = c(beta1, beta2)
sigma = 12

coefmat = NULL
infomat = NULL
for (sim in 1:50) {

	cat('sim:', sim, '\n')

	set.seed(sim*100+12345)

	n = 200

	## generate training data
	x1 = matrix(rnorm(n*q), ncol=q)  # generate a n-by-q matrix
	x2 = t(apply(x1, 1, combn, 2, prod)) # 1-row, 2-col, ????¨°?DD¡ê?¨¬?¨¢???¡Á?3??y¡ê?¨¢D3?¨¤¡ä¡ê?¡ã¡ä¨¢D??D¨°¨¢D3?¨¤¡ä
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)  #??? would a large sigma matter ???
	x1.tr = x1
	y.tr = y   # generate separate rather than together and then divide makes no difference

	## generate validation data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)
	x1.val = x1
	y.val = y

	## generate test data
	n = 10000   #??? the number of obs in test is large ???
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x = cbind(x1,x2)
	## sqrt(var(x%*%beta)/4)
	y = x%*%beta + rnorm(n,0,sigma)  # independent, so does not matter
	x1.te = x1
	y.te = y

	msevec = 0
	for (t in 1:9) {   # cross-validation for choosing "s"?
		
		fit = si(x=x1.tr, y=y.tr, s=0.1*t)  #??? what is s? Shrinkage ???
		
		## predict on validation set
		newx1.val = scale(x1.val, fit$mux, fit$sdx)  # use training mean and sd to scale validation
		newx2.val = t(apply(newx1.val, 1, combn, 2, prod))
		newx.val = cbind(newx1.val, newx2.val)
		newy.val = as.vector(y.val - fit$muy)  # why not y.val mean itself?
		yhat.val = newx.val%*%fit$bhat+fit$ahat

		msevec[t] = min(apply((yhat.val-newy.val)^2,2,mean))
	}
	t.opt = order(msevec)[1]  # order returns the index of ordering
	fit = si(x=x1.tr, y=y.tr, s=0.1*t.opt)
	yhat.val = newx.val%*%fit$bhat+fit$ahat
	idx = order(apply((yhat.val-newy.val)^2,2,mean))[1]
	ahat = fit$ahat[idx]
	bhat = fit$bhat[,idx]

	## calculate mse on test set
	newx1.te = scale(x1.te, fit$mux, fit$sdx)
	newx2.te = t(apply(newx1.te, 1, combn, 2, prod))
	newx.te = cbind(newx1.te, newx2.te)
	newy.te = as.vector(y.te - fit$muy)
	yhat.te = ahat + newx.te%*%bhat
	mse.te = mean((yhat.te-newy.te)^2)

	write.table(matrix(c(ahat,bhat),nrow=1),'coefmat.mycase2.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
	write.table(matrix(c(t.opt, fit$lambda[idx], mse.te),nrow=1),'infomat.mycase2.txt', col.names=FALSE, row.names=FALSE, append=TRUE)
}