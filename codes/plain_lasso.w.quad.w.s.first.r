## variable selection with interaction and quadratic terms
## formulation: OLS + lambda(s*main + (1-s)*interaction)
## the input contains all effects.
## update on 11/16/2017

plain.lasso.poly.sfirst = function(x, y, q, s, lambda=NULL, beta_ori=FALSE) {
	# user can specify lambda, otherwise glmnet default.
	# user can choose if original betas are returned.

	muy = mean(y)
	newy = y-muy
	
	## take care of two tuning parameters for main and interaction effects
	x1 = x[,1:q]/s
	x2 = x[,(q+1):(q*(q+3)/2)]/(1-s)
	x.temp = cbind(x1,x2)

	x.reg.std = scale(x.temp)
	mux = as.vector(attributes(x.reg.std)$'scaled:center')
	sdx = as.vector(attributes(x.reg.std)$'scaled:scale')
	# now all x's are standardized
	
	## fit lasso
	fit = glmnet(x.reg.std, newy, lambda=lambda, standardize=FALSE, intercept=TRUE)
	ahat = as.vector(fit$a0)
	bhat = as.matrix(fit$beta)

	## scale intercept and beta back to original
	if (beta_ori==TRUE) {
		bhat_ori = apply(bhat, 2, "/", sdx)
		ahat_ori = ahat - colSums(apply(bhat_ori, 2, "*", mux)) + muy
		bhat_ori[1:q,] = bhat_ori[1:q,]/s
		bhat_ori[-(1:q),] = bhat_ori[-(1:q),]/(1-s)
		return(list(ahat=ahat, bhat=bhat, ahat_ori = ahat_ori, bhat_ori = bhat_ori, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	} else {
		return(list(ahat=ahat, bhat=bhat, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	}
}

