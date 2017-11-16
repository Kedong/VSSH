rm(list=ls())

ptm <- proc.time()

setwd('D:/KC_file/Dropbox/VSSH')

library(glmnet)
library(ncvreg)

source('codes/vanilla_stepwise.w.quad.r')
#source('codes/vanilla_lasso.w.quad.r')
source('codes/vanilla_lasso.w.quad.double.std.r')
source('codes/shim.w.quad.kversion.r')
#source('codes/plain_lasso.w.quad.w.s.later.r')
source('codes/plain_lasso.w.quad.w.s.first.r')

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numintzero = 0  # number of interactions that were set as 0
numnonusemain = 0  # number of mains that dont have interactions
constcoef = 1  # if the coefficients are constant

main.magnitude = 4
int.magnitude = 1
quad.magnitude = 1

beta1 = c(rep(1,numsig),rep(0,q-numsig))
beta2 = combn(c(rep(1,numsig-numnonusemain),rep(0,q-numsig+numnonusemain)),2,prod)
if (numintzero>0) beta2[sample(which(beta2!=0),numintzero)]=0
beta3 = c(rep(1,numsig),rep(0,q-numsig))
beta1 = beta1 * main.magnitude
beta2 = beta2 * int.magnitude
beta3 = beta3 * quad.magnitude
beta = c(beta1, beta2, beta3)
if (constcoef!=1) {for (i in 1:length(beta)) beta[i]=beta[i]*(runif(1)+0.5)}
sigma = 8

sim.num = 50

#lambda.shim = exp(seq(20,-10,-1))
lambda.shim = exp(seq(15,-5,-0.5))
#lambda.shim = exp(seq(1,0,-0.5))  # this is for test
# let the algorithm choose lambda.lasso

coef.mat.lasso.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lasso.poly = numeric(sim.num)
coef.mat.lassoori.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lassoori.poly = numeric(sim.num)
coef.mat.step.poly = matrix(0,nrow=sim.num,ncol=q*(q+3)/2)
info.step.poly = numeric(sim.num)
coef.mat.stepori = matrix(0,nrow=sim.num,ncol=q*(q+3)/2)
info.stepori.poly = numeric(sim.num)
coef.mat.shim.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.shim.poly = numeric(sim.num)
signal2noise = numeric(sim.num)

for (sim in 1:sim.num) {

	cat('sim:', sim, '\n')

	set.seed(sim*100+12345)

	n = 200

	## generate training data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x3 = x1^2
	x.tr = cbind(x1,x2,x3)
	## sqrt(var(x.tr%*%beta)/4)
	y = x.tr%*%beta + rnorm(n,0,sigma)
	x1.tr = x1
	y.tr = y

	## generate validation data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x3 = x1^2
	x.val = cbind(x1,x2,x3)
	## sqrt(var(x.val%*%beta)/4)
	y = x.val%*%beta + rnorm(n,0,sigma)
	x1.val = x1
	y.val = y

	## generate test data
	n = 10000
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x3 = x1^2
	x.te = cbind(x1,x2,x3)
	## sqrt(var(x.te%*%beta)/4)
	y = x.te%*%beta + rnorm(n,0,sigma)
	x1.te = x1
	y.te = y

	signal2noise[sim] = var(x.te%*%beta)/sigma^2
	
	## SHIM Polynomial ##
#	coef.all.shim.poly.temp = matrix(,nrow=length(beta)+1,ncol=9)
#	for (t in 1:9) {
#		fit.shim.poly = shim.poly(x=x1.tr, y=y.tr, lambda = lambda.shim, t=t)
#		
#		## predict on validation set
#		yval.pred = matrix(rep(as.vector(fit.shim.poly$intercept),each=nrow(x.val)),ncol=length(lambda.shim)) +
#			x.val%*%fit.shim.poly$coefs
#		msevec.shim.poly = colMeans((apply(yval.pred,2,"-",y.val))^2)
#
#		coef.all.shim.poly.temp[,t] = c(fit.shim.poly$intercept[order(msevec.shim.poly)[1]],fit.shim.poly$coefs[,order(msevec.shim.poly)[1]])
#	}
#	yhat.val = cbind(1,x.val)%*%coef.all.shim.poly.temp
#	coef.mat.shim.poly[sim,] = coef.all.shim.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]][-1]
	## mse on test data ##
#	yhat.test = cbind(1,x.te)%*%matrix(coef.all.shim.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]],ncol=1)
#	info.shim.poly[sim] = mean((yhat.test-y.te)^2)
	
	
	## Our LASSO Polynomial ##
#	coef.all.lasso.poly.temp = matrix(,nrow=length(beta)+1,ncol=9)
#	for (t in 1:9) {
#		fit.lasso.poly = si.lasso.poly(x=x1.tr, y=y.tr, s=0.1*t, beta_ori=T)
		
		## predict on validation set
#		yval.pred = matrix(rep(fit.lasso.poly$ahat_ori,each=nrow(x.val)),ncol=length(fit.lasso.poly$ahat_ori)) +
#			x.val%*%fit.lasso.poly$bhat_ori
#		msevec.lasso.poly = colMeans((apply(yval.pred,2,"-",y.val))^2)
		
#		coef.all.lasso.poly.temp[,t] = c(fit.lasso.poly$ahat_ori[order(msevec.lasso.poly)[1]],fit.lasso.poly$bhat_ori[,order(msevec.lasso.poly)[1]])
#	}
#	yhat.val = cbind(1,x.val)%*%coef.all.lasso.poly.temp
#	coef.mat.lasso.poly[sim,] = coef.all.lasso.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]][-1]
	## mse on test data ##
#	yhat.test = cbind(1,x.te)%*%matrix(coef.all.lasso.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]],ncol=1)
#	info.lasso.poly[sim] = mean((yhat.test-y.te)^2)

	## Our LASSO Polynomial - double standardization ##

	fit.lasso.poly = si.lasso.poly.db(x=x1.tr, y=y.tr, beta_ori=T)
		
	## predict on validation set
	yval.pred = matrix(rep(fit.lasso.poly$ahat_ori,each=nrow(x.val)),ncol=length(fit.lasso.poly$ahat_ori)) +
			x.val%*%fit.lasso.poly$bhat_ori
	msevec.lasso.poly = colMeans((apply(yval.pred,2,"-",y.val))^2)
		
	coef.mat.lasso.poly[sim,] = fit.lasso.poly$bhat_ori[,order(msevec.lasso.poly)[1]]
	## mse on test data ##
	yhat.test = as.vector(x.te%*%fit.lasso.poly$bhat_ori[,order(msevec.lasso.poly)[1]]) + fit.lasso.poly$ahat_ori[order(msevec.lasso.poly)[1]]
	info.lasso.poly[sim] = mean((yhat.test-y.te)^2)
	
	## Our stepwise Polynomial ##
#	fit.ourstep.poly = si.stepwise.poly(x1.tr,y.tr)
#	coef.mat.step.poly[sim,] = fit.ourstep.poly$bhat_ori
	## mse on test data ##
#	yhat.test = cbind(1,x.te)%*%matrix(c(fit.ourstep.poly$ahat_ori,fit.ourstep.poly$bhat_ori),ncol=1)
#	info.step.poly[sim] = mean((yhat.test-y.te)^2)

	## LASSO 0.777##
	coef.all.plain.lasso.poly.temp = matrix(,nrow=length(beta)+1,ncol=9)
	for (t in 1:9) {
		fit.plain.lasso = plain.lasso.poly.sfirst(x=x.tr, y=y.tr, q=q, s=0.1*t, beta_ori=T)
		
		## predict on validation set
		yval.plain.pred = matrix(rep(fit.plain.lasso$ahat_ori,each=nrow(x.val)),ncol=length(fit.plain.lasso$ahat_ori)) +
			x.val%*%fit.plain.lasso$bhat_ori
		msevec.plain.lasso.poly = colMeans((apply(yval.plain.pred,2,"-",y.val))^2)
		
		coef.all.plain.lasso.poly.temp[,t] = c(fit.plain.lasso$ahat_ori[order(msevec.plain.lasso.poly)[1]],
						fit.plain.lasso$bhat_ori[,order(msevec.plain.lasso.poly)[1]])
	}

	yhat.plain.val = cbind(1,x.val)%*%coef.all.plain.lasso.poly.temp
	coef.mat.lassoori.poly[sim,] = coef.all.plain.lasso.poly.temp[,order(colMeans((apply(yhat.plain.val,2,"-",y.val))^2))[1]][-1]

	yhat.te = cbind(1,x.te)%*%matrix(coef.all.plain.lasso.poly.temp[,order(colMeans((apply(yhat.plain.val,2,"-",y.val))^2))[1]],ncol=1)
	info.lassoori.poly[sim] = mean((yhat.te-y.te)^2)  #mse for sim = 71.91585, t=5 output same as w/o s

	## 22222 LASSO 22222 0.816##

#	fit.plain.lasso = plain.lasso.poly.wos(x=x.tr, y=y.tr, q=q, beta_ori=T)
#	yval.plain.val = matrix(rep(fit.plain.lasso$ahat_ori,each=nrow(x.val)),ncol=length(fit.plain.lasso$ahat_ori)) +
#			x.val%*%fit.plain.lasso$bhat_ori
#	idx = order(colMeans(apply(yval.plain.val, 2, "-", y.val)^2))[1]
#	ahat.lasso = as.numeric(fit.plain.lasso$ahat_ori[idx])
#	bhat.lasso = as.vector(fit.plain.lasso$bhat_ori[,idx])
#	coef.mat.lassoori.poly[sim,] = matrix(bhat.lasso,nrow=1)

#	yhat.te = ahat.lasso + x.te%*%bhat.lasso
#	info.lassoori.poly[sim] = mean((yhat.te-y.te)^2)   #mse for sim1 = 72.33127

	## Stepwise (AIC) ##
	## Use only training to fit stepwise ##
	## Intercept included ##
#	data.step = cbind.data.frame(x.tr,y.tr)
#	tr.m0 = lm(y.tr~1, data=data.step)
#	tr.m1 = lm(y.tr~., data=data.step)
#	fit.step.poly = step(tr.m0, scope=list(upper=tr.m1, lower=tr.m0), direction="both", trace=0, data=data.step)
#	coef.mat.stepori[sim,as.numeric(names(attr(fit.step.poly$terms,"dataClasses"))[-1])] = as.numeric(fit.step.poly$coef[-1])
#	yhat.te.step = as.numeric(fit.step.poly$coef[1]) + x.te%*%coef.mat.stepori[sim,]
#	info.stepori.poly[sim] = mean((yhat.te.step-y.te)^2)
}


final.info = cbind(info.lasso.poly, info.lassoori.poly, info.step.poly, info.stepori.poly, info.shim.poly)
colnames(final.info) = c("mse.vanillalasso.poly", "mse.plainlasso.poly", "mse.vanillastep.poly", "mse.plainstep.poly", "mse.shim.poly")
final.info

snr.result = mean(signal2noise)

write.csv(t(coef.mat.shim.poly), "results/cond8/coef.plain.shim.csv")
write.csv(t(coef.mat.lasso.poly), "results/cond8/coef.vanilla.lasso.csv")
write.csv(t(coef.mat.lassoori.poly), "results/cond8/coef.plain.lasso.csv")
write.csv(t(coef.mat.step.poly), "results/cond8/coef.vanilla.stepwise.csv")
write.csv(t(coef.mat.stepori), "results/cond8/coef.plain.stepwise.csv")
write.csv(final.info, "results/cond8/mse.csv")
write.csv(snr.result, "results/cond8/snr.csv")

snr.result

proc.time() - ptm
