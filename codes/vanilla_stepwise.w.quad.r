## variable selection with interaction based upon stepwise
## Input only the main effects

si.stepwise.poly = function(x, y) {
	
	p = dim(x)[2]

	muy = mean(y)
	newy = y-muy

	x1 = scale(x)
	mux = as.vector(attributes(x1)$'scaled:center')
	sdx = as.vector(attributes(x1)$'scaled:scale')

	## create interaction term
	x2 = t(apply(x1, 1, combn, 2, prod))
	x3 = x1^2
	x.new = cbind(x1,x2,x3)
	newdata = cbind.data.frame(x.new,newy)

	## fit stepwise
	tr.m0 = lm(newy~1, data=newdata)
	tr.m1 = lm(newy~., data=newdata)
	fit.step = step(tr.m0, scope=list(upper=tr.m1, lower=tr.m0), direction="both", trace=0, data=newdata)
	bhat = rep(0,p*(p+1)/2+p)
	bhat[as.numeric(names(attr(fit.step$terms,"dataClasses"))[-1])] = as.numeric(fit.step$coef[-1])
	ahat = as.numeric(fit.step$coef[1])

	## scale intercept and beta back to original
	ahat_ori = ahat
	bhat_ori = bhat/c(sdx,combn(sdx,2,prod),sdx^2)
	bhat_final = bhat_ori

	bhat_int = numeric(p)
	for (i in 1:p){
		xmulti = rep(1,p)
		xmulti[i] = 0
		mux[i] = 1
		a = combn(mux,2,prod) * (rep(1,p*(p-1)/2) - combn(xmulti,2,prod))
		bhat_int[i] = sum(bhat_ori[(p+1):(p*(p+1)/2)]*a)
		mux = as.vector(attributes(x1)$'scaled:center')
	}
	bhat_final[1:p] = bhat_ori[1:p] - bhat_int - 2*(bhat_ori[(p*(p+1)/2+1):(p*(p+1)/2+p)]*mux)
	ahat_final = ahat_ori - sum(bhat_ori[1:p]*mux) + sum(bhat_ori[(p+1):(p*(p+1)/2)]*combn(mux,2,prod)) +
		+ sum(bhat_ori[(p*(p+1)/2+1):(p*(p+1)/2+p)]*mux^2) + muy

	return(list(ahat = ahat, bhat = bhat, ahat_ori = ahat_final, bhat_ori = bhat_final, mux=mux, sdx=sdx, muy=muy))
}

