## variable selection with interaction and quadratic terms
## formulation: OLS + lambda(s*main + (1-s)*interaction)
## the input contains all effects.
## update on 11/16/2017

plain.lasso.poly = function(x, y, q, s, lambda=NULL, beta_ori=FALSE) {
	# user can specify lambda, otherwise glmnet default.
	# user can choose if original betas are returned.

	muy = mean(y)
	newy = y-muy

	x.reg.std = scale(x)
	mux = as.vector(attributes(x.reg.std)$'scaled:center')
	sdx = as.vector(attributes(x.reg.std)$'scaled:scale')
	# now all x's are standardized
	
	## take care of two tuning parameters for main and interaction effects
	x1 = x.reg.std[,1:q]/s
	x2 = x.reg.std[,(q+1):(q*(q+1)/2)]/(1-s)
	x3 = x.reg.std[,(q*(q+1)/2+1):(q*(q+3)/2)]/(1-s)

	newx = cbind(x1,x2,x3)

	## fit lasso
	fit = glmnet(newx, newy, lambda=lambda, standardize=FALSE, intercept=TRUE)
	ahat = as.vector(fit$a0)
	bhat = as.matrix(fit$beta)
	bhat[1:q,] = bhat[1:q,]/s
	bhat[-(1:q),] = bhat[-(1:q),]/(1-s)

	## scale intercept and beta back to original
	if (beta_ori==TRUE) {
		bhat_ori = apply(bhat, 2, "/", sdx)
		ahat_ori = ahat - colSums(apply(bhat_ori, 2, "*", mux)) + muy
		return(list(ahat=ahat, bhat=bhat, ahat_ori = ahat_ori, bhat_ori = bhat_ori, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	} else {
		return(list(ahat=ahat, bhat=bhat, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	}
}

