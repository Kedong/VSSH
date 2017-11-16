## variable selection with interaction
## formulation: OLS + lambda(s*main + (1-s)*interaction)
## the input is main effects only. the code creates all two-way interactions by itself.
## update on 12/20/2015

si.lasso.poly.db = function(x, y, lambda=NULL, beta_ori=FALSE) {
	# user can specify lambda, otherwise glmnet default.
	# user can choose if original betas are returned.
	
	p = dim(x)[2]

	muy = mean(y)
	newy = y-muy

	x1 = scale(x)
	mux = as.vector(attributes(x1)$'scaled:center')
	sdx = as.vector(attributes(x1)$'scaled:scale')

	## create interaction term
	x2 = t(apply(x1, 1, combn, 2, prod))
	x3 = x1^2

	newx = cbind(x1,x2,x3)

	## fit lasso
	fit = glmnet(newx, newy, lambda=lambda)
	ahat = as.vector(fit$a0)
	bhat = as.matrix(fit$beta)

	## scale intercept and beta back to original
	if (beta_ori==TRUE) {
		ahat_ori = ahat
		bhat_ori = apply(bhat,2,"/",c(sdx,combn(sdx,2,prod),sdx^2))
		bhat_final = bhat_ori

		bhat_int = matrix(,p,ncol(bhat))
		for (i in 1:p){
			xmulti = rep(1,p)
			xmulti[i] = 0
			mux[i] = 1
			a = combn(mux,2,prod) * (rep(1,p*(p-1)/2) - combn(xmulti,2,prod))
			bhat_int[i,] = colSums(apply(bhat_ori[(p+1):(p*(p+1)/2),],2,"*",a))
			mux = as.vector(attributes(x1)$'scaled:center')
		}
		bhat_final[1:p,] = bhat_ori[1:p,] - bhat_int - 2*apply(bhat_ori[(p*(p+1)/2+1):(p*(p+1)/2+p),],2,"*",mux)
		ahat_final = ahat_ori - as.vector(colSums(apply(bhat_ori[1:p,],2,"*",mux))) +
				as.vector(colSums(apply(bhat_ori[(p+1):(p*(p+1)/2),],2,"*",combn(mux,2,prod)))) +
				as.vector(colSums(apply(bhat_ori[(p*(p+1)/2+1):(p*(p+1)/2+p),],2,"*",mux^2))) + muy

		return(list(ahat=ahat, bhat=bhat, ahat_ori = ahat_final, bhat_ori = bhat_final, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	} else {
		return(list(ahat=ahat, bhat=bhat, mux=mux, sdx=sdx, muy=muy, lambda=fit$lambda))
	}
}

