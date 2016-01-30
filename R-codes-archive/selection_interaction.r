## variable selection with interaction
## formulation: OLS + lambda(s*main + (1-s)*interaction)
## the input is main effects only. the code creates all two-way interactions by itself.
## update on 10/26/2015

si = function(x, y, s, lambda) {
	
	n = dim(x)[1]
	p = dim(x)[2]

	muy = mean(y)
	newy = y-muy

	x1 = scale(x)
	mux = as.vector(attributes(x1)$`scaled:center`)
	sdx = as.vector(attributes(x1)$`scaled:scale`)
	x1 = x1

	## create interaction term
	x2 = t(apply(x1, 1, combn, 2, prod))
	
	## take care of two tuning parameters for main and interaction effects
	x1 = x1/s
	x2 = x2/(1-s)

	newx = cbind(x1,x2)

	## fit lasso
	fit = glmnet(newx, newy, standardize=FALSE, intercept=TRUE)
	ahat = as.vector(fit$a0)
	bhat = as.matrix(fit$beta)
	bhat[1:p,] = bhat[1:p,]/s
	bhat[-(1:p),] = bhat[-(1:p),]/(1-s)

	return(list(ahat=ahat, bhat=bhat, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
}