rm(list=ls())

ptm <- proc.time()

setwd('D:/KC_file/Dropbox/VSSH')

library(glmnet)
library(ncvreg)

source('codes/vanilla_lasso.w.quad.r')

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numintzero = 0  # number of interactions that were set as 0
numnonusemain = 0  # number of mains that dont have interactions
constcoef = 1  # if the coefficients are constant

main.magnitude = 5
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
sigma = 3

## part of robustness##
n.main.zero = 2
beta[1:n.main.zero]=0

sim.num = 50

coef.mat.lasso.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lasso.poly = numeric(sim.num)
coef.mat.lassoori.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lassoori.poly = numeric(sim.num)
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
	
	## Our LASSO Polynomial ##
	coef.all.lasso.poly.temp = matrix(,nrow=length(beta)+1,ncol=9)
	for (t in 1:9) {
		fit.lasso.poly = si.lasso.poly(x=x1.tr, y=y.tr, s=0.1*t, beta_ori=T)
		
		## predict on validation set
		yval.pred = matrix(rep(fit.lasso.poly$ahat_ori,each=nrow(x.val)),ncol=length(fit.lasso.poly$ahat_ori)) +
			x.val%*%fit.lasso.poly$bhat_ori
		msevec.lasso.poly = colMeans((apply(yval.pred,2,"-",y.val))^2)
		
		coef.all.lasso.poly.temp[,t] = c(fit.lasso.poly$ahat_ori[order(msevec.lasso.poly)[1]],fit.lasso.poly$bhat_ori[,order(msevec.lasso.poly)[1]])
	}
	yhat.val = cbind(1,x.val)%*%coef.all.lasso.poly.temp
	coef.mat.lasso.poly[sim,] = coef.all.lasso.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]][-1]
	## mse on test data ##
	yhat.test = cbind(1,x.te)%*%matrix(coef.all.lasso.poly.temp[,order(colMeans((apply(yhat.val,2,"-",y.val))^2))[1]],ncol=1)
	info.lasso.poly[sim] = mean((yhat.test-y.te)^2)
	

	## LASSO ##
	fit.lasso = glmnet(x=x.tr, y=y.tr)
	yhat.val = predict(fit.lasso, x.val)
	idx = order(apply((yhat.val-as.vector(y.val))^2,2,mean))[1]

	ahat.lasso = as.numeric(fit.lasso$a0[idx])
	bhat.lasso = as.vector(fit.lasso$beta[,idx])
	coef.mat.lassoori.poly[sim,] = matrix(bhat.lasso,nrow=1)

	yhat.te = ahat.lasso + x.te%*%bhat.lasso
	info.lassoori.poly[sim] = mean((yhat.te-y.te)^2)
}

snr.rob = mean(signal2noise)


final.info = cbind(info.lasso.poly, info.lassoori.poly)
colnames(final.info) = c("mse.vanillalasso.poly", "mse.plainlasso.poly")
final.info

write.csv(t(coef.mat.lasso.poly), "results/cond1.robust/coef.vanilla.lasso.csv")
write.csv(t(coef.mat.lassoori.poly), "results/cond1.robust/coef.plain.lasso.csv")
write.csv(final.info, "results/cond1.robust/mse.csv")
write.csv(snr.rob, "results/cond1.robust/snr.rob.csv")

proc.time() - ptm
