setwd('D:/Schooling/Higher Education/!Research/SHIM_New/VSSH')

library(glmnet)
library(ncvreg)

source('R-codes-use/selection_interaction_lasso.poly.r')
source('R-codes-use/proportion_heredity.r')

gene = read.table("real data/gene_v1.txt",header=T)
# this version is to select people first, then select genes #
# see if the results are consistent with previous #
gene = read.table("real data/gene_v2.txt",header=T)
# this version is to select genes first, then select people #
gene_pure = gene[,-(1:2)]
day2death = log(gene[,2])
n_gene = 100
q = n_gene
total_beta = n_gene*2 + n_gene*(n_gene-1)/2

slt_idx = order(apply(gene_pure, 2, sd),decreasing=T)
gene_x.1 = gene_pure[,slt_idx[1:n_gene]]
gene_x.2 = t(apply(gene_x.1, 1, combn, 2, prod))
gene_x.3 = gene_x.1^2
gene_x = as.matrix(cbind(gene_x.1,gene_x.2,gene_x.3))

sim.num = 50
tr_portion = 0.6
val_portion = 0.2
te_portion = 0.2

coef.mat.lasso.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lasso.poly = numeric(sim.num)
coef.mat.lassoori.poly = matrix(,nrow=sim.num,ncol=q*(q+3)/2)
info.lassoori.poly = numeric(sim.num)
signal2noise = numeric(sim.num)

for (sim in 1:sim.num) {

	cat('sim:', sim, '\n')
	set.seed(sim*100+12345)

	# use 60%, 20%, 20% to generate data #
	tr_idx = sample(nrow(gene_x),round(nrow(gene_x)*tr_portion))
	x1.tr = gene_x[tr_idx,1:n_gene]
	x.tr = gene_x[tr_idx,]
	y.tr = day2death[tr_idx]
	x.temp = gene_x[-tr_idx,]
	y.temp = day2death[-tr_idx]

	val_idx = sample(nrow(x.temp),round(nrow(x.temp)*(val_portion/(1-tr_portion))))	
	x1.val = x.temp[val_idx,1:n_gene]
	x.val = x.temp[val_idx,]
	y.val = y.temp[val_idx]
	x1.te = x.temp[-val_idx,1:n_gene]
	x.te = x.temp[-val_idx,]
	y.te = y.temp[-val_idx]

	## Our LASSO Polynomial ##
	coef.all.lasso.poly.temp = matrix(,nrow=total_beta+1,ncol=9)
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

final.info = cbind(info.lasso.poly, info.lassoori.poly)
colnames(final.info) = c("mse.olasso.poly", "mse.lasso.ori.poly")
final.info

