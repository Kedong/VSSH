rm(list=ls())

ptm <- proc.time()

setwd('D:/Schooling/Higher Education/!papers&research/SHIM_New')
library(glmnet)
source('selection_interaction_w-beta_ori.r')
source('shim.r')

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numintzero = 0  # number of interactions that were set as 0
numnonusemain = 0  # number of mains that dont have interactions
constcoef = 1  # if the coefficients are constant
main.magnitude = 0.1
int.magnitude = 10

beta1 = c(rep(1,numsig),rep(0,q-numsig))
beta2 = combn(c(rep(1,numsig-numnonusemain),rep(0,q-numsig+numnonusemain)),2,prod)
if (numintzero>0) beta2[sample(which(beta2!=0),numintzero)]=0
beta1 = beta1 * main.magnitude
beta2 = beta2 * int.magnitude
beta = c(beta1, beta2)
if (constcoef!=1) {for (i in 1:length(beta)) beta[i]=beta[i]*(runif(1)+0.5)}
sigma = 12

sim.num = 6
coef.mat.our = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.our = matrix(,sim.num,3)
coef.mat.lasso = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.lasso = matrix(,sim.num,1)
coef.mat.step = matrix(0,nrow=sim.num,ncol=q*(q+1)/2)
info.step = matrix(,sim.num,1)
coef.mat.shim = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.shim = matrix(,sim.num,1)

for (sim in 1:sim.num) {

	cat('sim:', sim, '\n')

	set.seed(sim*100+12345)

	n = 200

	## generate training data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x.tr = cbind(x1,x2)
	## sqrt(var(x.tr%*%beta)/4)
	y = x.tr%*%beta + rnorm(n,0,sigma)
	x1.tr = x1
	y.tr = y

	## generate validation data
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x.val = cbind(x1,x2)
	## sqrt(var(x.val%*%beta)/4)
	y = x.val%*%beta + rnorm(n,0,sigma)
	x1.val = x1
	y.val = y

	## generate test data
	n = 10000
	x1 = matrix(rnorm(n*q), ncol=q)
	x2 = t(apply(x1, 1, combn, 2, prod))
	x.te = cbind(x1,x2)
	## sqrt(var(x.te%*%beta)/4)
	y = x.te%*%beta + rnorm(n,0,sigma)
	x1.te = x1
	y.te = y

	## Non-adaptive SHIM ##
	## Transform lambda1&2 to lambda & s ##
	## phi here is the coefs with intercept ##
	lambda.list = c(1, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 125, 150, 175, 200, 225, 250) #185 does not converge...
	lambda.list.length = length(lambda.list)
	msevec.shim.1 = 0
	v.list = 0
	for (k in 1:lambda.list.length) {
		msevec.shim.2 = 0
		for (v in 1:9) {
			fit.shim = shim(x=x1.tr, y=y.tr, lambda1=lambda.list[k]*v*0.1, lambda2=lambda.list[k]*(1-v*0.1))
			## predict on validation set
			yval.pred = fit.shim$phi[1] + x.val%*%fit.shim$phi[-1]
			msevec.shim.2[v] = mean((yval.pred-y.val)^2)
		}
		v.opt = order(msevec.shim.2)[1]
		fit.shim = shim(x=x1.tr, y=y.tr, lambda1=lambda.list[k]*v.opt*0.1, lambda2=lambda.list[k]*(1-v.opt*0.1))
		yval.pred = fit.shim$phi[1] + x.val%*%fit.shim$phi[-1]
		msevec.shim.1[k] = mean((yval.pred-y.val)^2)
		v.list[k] = v.opt
	}
	k.opt = order(msevec.shim.1)[1]
	fit.shim = shim(x=x1.tr, y=y.tr, lambda1=lambda.list[k.opt]*v.list[k.opt]*0.1, lambda2=lambda.list[k.opt]*(1-v.list[k.opt]*0.1))
	coef.mat.shim[sim,] = fit.shim$phi[-1]
	## calculate mse on test set
	yte.pred = fit.shim$phi[1] + x.te%*%fit.shim$phi[-1]
	info.shim[sim,] = mean((yte.pred-y.te)^2)

}

colnames(info.shim) = "mse.shim"

## Sensitivity and Specificity ##

cnt.criterion = main.magnitude/100
#cnt.criterion = self-determine
se.our = apply(abs(coef.mat.our[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.our = apply(abs(coef.mat.our[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.lasso = apply(abs(coef.mat.lasso[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.lasso = apply(abs(coef.mat.lasso[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.step = apply(abs(coef.mat.step[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.step = apply(abs(coef.mat.step[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.shim = apply(abs(coef.mat.shim[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.shim = apply(abs(coef.mat.shim[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)

se.our.avg = mean(se.our)
se.our.sd = sd(se.our)
sp.our.avg = mean(sp.our)
sp.our.sd = sd(sp.our)
se.lasso.avg = mean(se.lasso)
se.lasso.sd = sd(se.lasso)
sp.lasso.avg = mean(sp.lasso)
sp.lasso.sd = sd(sp.lasso)
se.step.avg = mean(se.step)
se.step.sd = sd(se.step)
sp.step.avg = mean(sp.step)
sp.step.sd = sd(sp.step)
se.shim.avg = mean(se.shim)
se.shim.sd = sd(se.shim)
sp.shim.avg = mean(sp.shim)
sp.shim.sd = sd(sp.shim)

final.se.sp = cbind(se.our.avg,sp.our.avg,se.lasso.avg,sp.lasso.avg,se.step.avg,sp.step.avg,se.shim.avg,sp.shim.avg,
	se.our.sd,sp.our.sd,se.lasso.sd,sp.lasso.sd,se.step.sd,sp.step.sd,se.shim.sd,sp.shim.sd)
final.se.sp

#write.csv(final.info,"mse.comparison.cond1.csv",row.names=F)
#write.table(t(final.se.sp),"final.cond1.csv",col.names=F,sep=",")

proc.time() - ptm
