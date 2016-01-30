rm(list=ls())

ptm <- proc.time()

setwd('D:/Schooling/Higher Education/!papers&research/SHIM_New')
library(glmnet)
source('selection_interaction_w-beta_ori.r')
source('shim.r')
source('selection_interaction_stepwise.r')

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numintzero = 0  # number of interactions that were set as 0
numnonusemain = 0  # number of mains that dont have interactions
constcoef = 1  # if the coefficients are constant
main.magnitude = 1
int.magnitude = 1

beta1 = c(rep(1,numsig),rep(0,q-numsig))
beta2 = combn(c(rep(1,numsig-numnonusemain),rep(0,q-numsig+numnonusemain)),2,prod)
if (numintzero>0) beta2[sample(which(beta2!=0),numintzero)]=0
beta1 = beta1 * main.magnitude
beta2 = beta2 * int.magnitude
beta = c(beta1, beta2)
if (constcoef!=1) {for (i in 1:length(beta)) beta[i]=beta[i]*(runif(1)+0.5)}
sigma = 8

sim.num = 100
coef.mat.our = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.our = matrix(,sim.num,3)
coef.mat.lasso = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.lasso = matrix(,sim.num,1)
coef.mat.step = matrix(0,nrow=sim.num,ncol=q*(q+1)/2)
info.step = matrix(,sim.num,1)
coef.mat.ourstep = matrix(,nrow=sim.num,ncol=q*(q+1)/2)
info.ourstep = matrix(,sim.num,1)
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
	
	## OUR METHOD - LASSO ##
	msevec.our = 0
	#for (t in 1:9) {	
	#	fit = si(x=x1.tr, y=y.tr, s=0.1*t, beta_ori=T)
	#	
	#	## predict on validation set
	#	y_ori_fit = x.val %*% fit$bhat_ori + as.vector(fit$ahat_ori)
	#	msevec[t] = min(apply((y_ori_fit-as.vector(y.val))^2,2,mean))
	#}

	# MSE slightly differs. The order is the same.
	for (t in 1:9) {
		
		fit.our = si(x=x1.tr, y=y.tr, s=0.1*t, beta_ori=F)
		
		## predict on validation set
		newx1.val = scale(x1.val, fit.our$mux, fit.our$sdx)
		newx2.val = t(apply(newx1.val, 1, combn, 2, prod))
		newx.val = cbind(newx1.val, newx2.val)
		newy.val = as.vector(y.val - fit.our$muy)
		yhat.val = newx.val%*%fit.our$bhat+as.vector(fit.our$ahat)
		
		msevec.our[t] = min(apply((yhat.val-newy.val)^2,2,mean))
	}
	t.opt = order(msevec.our)[1]
	fit.our = si(x=x1.tr, y=y.tr, s=0.1*t.opt, beta_ori=T)
	yhat.val = newx.val%*%fit.our$bhat+fit.our$ahat
	idx = order(apply((yhat.val-newy.val)^2,2,mean))[1]
	ahat.our = fit.our$ahat_ori[idx]
	bhat.our = fit.our$bhat_ori[,idx]

	## calculate mse on test set
	newx1.te = scale(x1.te, fit.our$mux, fit.our$sdx)
	newx2.te = t(apply(newx1.te, 1, combn, 2, prod))
	newx.te = cbind(newx1.te, newx2.te)
	newy.te = as.vector(y.te - fit.our$muy)
	yhat.te = fit.our$ahat[idx] + newx.te%*%fit.our$bhat[,idx]
	mse.te.our = mean((yhat.te-newy.te)^2)
	
	coef.mat.our[sim,] = matrix(bhat.our,nrow=1)
	info.our[sim,] = matrix(c(t.opt, fit.our$lambda[idx], mse.te.our),nrow=1)

	## LASSO ##
	fit.lasso = glmnet(x=x.tr, y=y.tr)
	yhat.val = predict(fit.lasso, x.val)
	idx = order(apply((yhat.val-as.vector(y.val))^2,2,mean))[1]

	ahat.lasso = as.numeric(fit.lasso$a0[idx])
	bhat.lasso = as.vector(fit.lasso$beta[,idx])
	coef.mat.lasso[sim,] = matrix(bhat.lasso,nrow=1)

	yhat.te = ahat.lasso + x.te%*%bhat.lasso
	info.lasso[sim,] = mean((yhat.te-y.te)^2)

	## Choose a model by AIC in a Stepwise Algorithm ##
	## Use only training to fit stepwise ##
	## Intercept included ##
	data.step = cbind.data.frame(x.tr,y.tr)
	tr.m0 = lm(y.tr~1, data=data.step)
	tr.m1 = lm(y.tr~., data=data.step)
	fit.step = step(tr.m0, scope=list(upper=tr.m1, lower=tr.m0), direction="both", trace=0, data=data.step)
	coef.mat.step[sim,as.numeric(names(attr(fit.step$terms,"dataClasses"))[-1])] = as.numeric(fit.step$coef[-1])
	yhat.te.step = as.numeric(fit.step$coef[1]) + x.te%*%coef.mat.step[sim,]
	info.step[sim,] = mean((yhat.te.step-y.te)^2)

	## Our method applied to stepwise ##
	fit.ourstep = si.stepwise(x1.tr,y.tr)
	coef.mat.ourstep[sim,] = fit.ourstep$bhat_ori
	yhat.te.ourstep = fit.ourstep$ahat_ori + x.te%*%fit.ourstep$bhat_ori
	info.ourstep[sim,] = mean((yhat.te.ourstep-y.te)^2)

}

colnames(info.our) = c("t.opt","lambda.opt","mse.our")
colnames(info.lasso) = "mse.lasso"
colnames(info.step) = "mse.step"
colnames(info.ourstep) = "mse.ourstep"
final.info = cbind(info.our,info.lasso,info.step,info.ourstep)

## Sensitivity and Specificity ##

cnt.criterion = main.magnitude/100
#cnt.criterion = self-determine
se.our = apply(abs(coef.mat.our[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.our = apply(abs(coef.mat.our[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.lasso = apply(abs(coef.mat.lasso[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.lasso = apply(abs(coef.mat.lasso[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.step = apply(abs(coef.mat.step[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.step = apply(abs(coef.mat.step[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)
se.ourstep = apply(abs(coef.mat.ourstep[,which(beta!=0)])>cnt.criterion,1,sum) /sum(beta!=0)
sp.ourstep = apply(abs(coef.mat.ourstep[,which(beta==0)])<=cnt.criterion,1,sum) /sum(beta==0)

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
se.ourstep.avg = mean(se.ourstep)
se.ourstep.sd = sd(se.ourstep)
sp.ourstep.avg = mean(sp.ourstep)
sp.ourstep.sd = sd(sp.ourstep)

final.se.sp = cbind(se.our.avg,sp.our.avg,se.lasso.avg,sp.lasso.avg,se.step.avg,sp.step.avg,se.ourstep.avg,sp.ourstep.avg,
	se.our.sd,sp.our.sd,se.lasso.sd,sp.lasso.sd,se.step.sd,sp.step.sd,se.ourstep.sd,sp.ourstep.sd)
final.se.sp

write.csv(final.info,"mse.comparison.cond3.csv",row.names=F)
write.table(t(final.se.sp),"final.cond3.csv",col.names=F,sep=",")

proc.time() - ptm
