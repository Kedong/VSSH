########## SHIM for polynomial ##########
##### Kedong Chen, chen3548@umn.edu #####

shim.poly <- function(x, y, lambda, t, coef.prev=NULL, max.iter=10^4, eps1 = 10^(-10), eps2 = 10^(-4)){
	
	require(glmnet)

	np <- dim(x)
	n <- np[1]
	p <- np[2]
	lambda.beta = lambda*t*0.1
	lambda.gamma = lambda*(1-t*0.1)   # lambda is a vector

	## calculating interaction & poly terms ##
	xx = t(apply(x, 1, combn, 2, prod))
	xsq = x^2
	y.centered = y - mean(y)
	x.norm = scale(x)
	xx.norm = scale(xx)
	xsq.norm = scale(xsq)


	## Initial Estimate ##
	## coef.prev needs to be #predictor*#lambda ##
	if(is.null(coef.prev)){
		m0 = lm(y.centered~cbind(x.norm,xx.norm,xsq.norm))
		coef.x = m0$coef[2:(p+1)]
		coef.x[ abs(coef.x) <= eps1 ] <- eps1
		coef.xx = m0$coef[(p+2):(p*(p-1)/2+p+1)]
		coef.xsq = m0$coef[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1)]
		
		coef.prev = matrix(rep(as.vector(c(m0$coef[1],coef.x,coef.xx,coef.xsq)),length(lambda)),ncol=length(lambda))
	}

	coef.new = matrix(,nrow=nrow(coef.prev),ncol=ncol(coef.prev))
	iter.result = numeric(length(lambda))

	##### Iterative Approximation ######
	for (lam.ind in 1:length(lambda)){
	
	iter <- 0
	differ <- 1
	crit.old <- crossprod(y.centered - cbind(x.norm,xx.norm,xsq.norm)%*%coef.prev[-1,lam.ind]) +
		lambda.beta[lam.ind]*sum(abs(coef.prev[2:(p+1),lam.ind])) +
		lambda.gamma[lam.ind]*(sum(abs(coef.prev[(p+2):(p*(p-1)/2+p+1),lam.ind] / combn(coef.prev[2:(p+1),lam.ind], 2, prod))) +
			sum(abs(coef.prev[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1),lam.ind] / coef.prev[2:(p+1),lam.ind])))

	while(differ > eps2){

    ################################
    ### STEP 1: solve for \gamma ###
   	################################
		y.new = y.centered - x.norm%*%coef.prev[2:(p+1),lam.ind]
		xx.new = t(apply(xx.norm,1,"*",combn(coef.prev[2:(p+1),lam.ind], 2, prod)))
		xsq.new = t(apply(xsq.norm,1,"*",as.vector(coef.prev[2:(p+1),lam.ind])))
		#gamma = glmnet(x=cbind(xx.new,xsq.new), y=y.new, lambda = lambda.gamma[lam.ind], standardize=F)$beta
		gamobj = solve.gam(cbind(xx.new,xsq.new), y.new, lambda.gamma[lam.ind], rep(1,(p+1)*p/2))
		gamma = gamobj$phi.new

  	###############################
  	### STEP 2: solve for \beta ###
  	###############################
		for (j in 1:p){
			xmulti = rep(1,p)
			xmulti[j] = 0
			int.index = combn(xmulti,2,prod)

			y.new = y.centered - x.norm[,-j]%*%coef.prev[2:(p+1),lam.ind][-j] -
				as.matrix(xx.norm[,which(int.index!=0)]) %*% 
				(combn(coef.prev[2:(p+1),lam.ind], 2, prod)[which(int.index!=0)]*as.vector(gamma)[which(int.index!=0)]) -
				(xsq.norm[,-j]) %*% ((coef.prev[2:(p+1),lam.ind][-j]) * (as.vector(gamma)[(p*(p-1)/2+1):(p*(p-1)/2+p)][-j]))
			
			x.new = x.norm[,j] +
				as.matrix(xx.norm[,which(int.index==0)]) %*%
				((combn(coef.prev[2:(p+1),lam.ind], 2, prod)[which(int.index==0)]/(coef.prev[2:(p+1),lam.ind][j]))*as.vector(gamma)[which(int.index==0)]) +
				as.vector(gamma)[(p*(p-1)/2+1):(p*(p-1)/2+p)][j]*xsq.norm[,j]

			## coordinate descent ##
			b.hat = lm(y.new~x.new)$coef[2]
			if (lambda.beta[lam.ind]/2 < crossprod(x.new,y.new)){
				coef.new[j+1,lam.ind] = (crossprod(x.new,y.new)-lambda.beta[lam.ind]/2)/crossprod(x.new)
			} else if (lambda.beta[lam.ind]/(-2) > crossprod(x.new,y.new)){
				coef.new[j+1,lam.ind] = (crossprod(x.new,y.new)+lambda.beta[lam.ind]/2)/crossprod(x.new)
			} else {
				coef.new[j+1,lam.ind] = eps1
			}
		}
		coef.new[(p+2):(p*(p-1)/2+p+1),lam.ind] = combn(coef.new[2:(p+1),lam.ind], 2, prod) * as.vector(gamma)[1:(p*(p-1)/2)]
		coef.new[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1),lam.ind] = coef.new[2:(p+1),lam.ind] * as.vector(gamma)[(p*(p-1)/2+1):(p*(p-1)/2+p)]

		crit.new <- crossprod(y.centered - cbind(x.norm,xx.norm,xsq.norm)%*%coef.new[-1,lam.ind]) +
			lambda.beta[lam.ind]*sum(abs(coef.new[2:(p+1),lam.ind])) +
			lambda.gamma[lam.ind]*(sum(abs(as.vector(gamma))))

		differ <- abs(crit.new - crit.old) / abs(crit.old)
		crit.old <- crit.new
		coef.prev[,lam.ind] = coef.new[,lam.ind]
		iter <- iter + 1
		if(iter >= max.iter){
			cat("Warning: the whole iteration did not converge.. \n")
			break
		}
	} ## end while
	iter.result[lam.ind] = iter
	## scale betas back ##
	coef.new[2:(p+1),lam.ind] = coef.new[2:(p+1),lam.ind]/attr(scale(x),"scaled:scale")
	coef.new[(p+2):(p*(p-1)/2+p+1),lam.ind] = coef.new[(p+2):(p*(p-1)/2+p+1),lam.ind]/attr(scale(xx),"scaled:scale")
	coef.new[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1),lam.ind] = coef.new[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1),lam.ind]/attr(scale(xsq),"scaled:scale")
	coef.new[1,lam.ind] = mean(y) - sum(coef.new[2:(p+1),lam.ind]*attr(scale(x),"scaled:center")) -
		sum(coef.new[(p+2):(p*(p-1)/2+p+1),lam.ind]*attr(scale(xx),"scaled:center"))-
		sum(coef.new[(p*(p-1)/2+p+2):(p*(p-1)/2+2*p+1),lam.ind]*attr(scale(xsq),"scaled:center"))

	} ## end for

	coef.new[abs(coef.new) <= eps1] = 0
	return( list(intercept=coef.new[1,], coefs=coef.new[-1,], iter=iter.result) )
}

## sub-funtion adpated from original SHIM code ##
solve.gam <- function(x, y, lambda, gam, max.iter=10^4, eps1=10^(-10), eps2=10^(-4)) {

  xtx <- t(x) %*% x
  xty <- t(x) %*% y

  ## Initial estimate ##
  phi.old <- xty / apply(x^2, 2, sum)

  Dphi <- diag( length(gam) )
  Qold <- crossprod(y - x%*%phi.old) + lambda*sum(abs(phi.old)/gam)

  differ <- 1
  iter <- 1
  #sing <- 0

  while (differ > eps2) {
    denom <- phi.old*gam
    denom[ abs(denom) < eps1 ] <- eps1
    diag(Dphi) <- as.vector( 1 / abs(denom) )
    xtx2 <- xtx + lambda * Dphi
    phi.new <- qr.solve(xtx2, tol=1e-10) %*% xty

    Qnew <- crossprod(y - x%*%phi.new) + lambda*sum(abs(phi.new)/gam)
    differ <- sum(abs(Qnew - Qold)) / sum(abs(Qold))

    Qold <- Qnew
    phi.old <- phi.new

    iter <- iter + 1

    if(iter >= max.iter){
      cat("... gamma did not converge ...\n")
      break
    }

  }
  return(list(phi.new=phi.new))
}