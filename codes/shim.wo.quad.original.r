##
## NOTE:
## Please see "Example" at the bottom before running the code.
##
##
## main function: shim
##
## ---- Input Variables ----
## x: a matrix containing main terms
## lambda1: regularization parameter for betas
## lambda2: regularization parameter for gammas
## phi.weight: a vector containing adaptive weights
## phi.prev: can plug in coefficient estimates obtained before if any.
##           e.g. when fitting a model on a grid of lambdas, 
##                can use the coeff estimates obtained for the previous pair of (lambda1, lambda2) for faster convergence.
##           This will be used as initial values for the iteration.
##           Default is to use the OLS estimates from simple linear regression models with each individual predictor.
##
## ---- Stopping criterion ----
## The iteration stops when
##   1) differ = abs(crit.new - crit.old) / abs(crit.old) < eps2
##      (when it converged)
##   OR,
##   2) the number of iterations reached at max.iter
##      (when it didn't converge by then.. an warning will appear in this case)
## 






##
## main function: shim
##

shim <- function(x, y, lambda1, lambda2, phi.weight=NULL, phi.prev=NULL, max.iter=10^4, eps1 = 10^(-10), eps2 = 10^(-4)){

  np <- dim(x)
  n <- np[1]
  p <- np[2]
  K <- (p-1)*p/2


  ## index set for interactions ##
  jk1 <- jk2 <- numeric(K)
  i <- 1
  for(j in 1:(p-1)){
    for(k in (j+1):p){
      jk1[i] <- j; jk2[i] <- k
      i <- i + 1
    }
  }

  ## calculating interaction terms ##
  xx = NULL
  for(j in 1:K){
    xx = cbind(xx, x[,jk1[j]]*x[,jk2[j]])
  }


  ## Center y ##
  meany <- mean(y)
  y <- y - meany

  ## Standardize x ##
  meanx <- apply(x, 2, mean)
  x <- scale(x, center=meanx, scale=F)
  normx <- sqrt(apply(x^2, 2, mean))
  x <- scale(x, center=F, scale=normx)

  meanxx <- apply(xx, 2, mean)
  xx <- scale(xx, center=meanxx, scale=F)
  normxx <- sqrt(apply(xx^2, 2, mean))
  xx <- scale(xx, center=F, scale=normxx)

  ## Initial Estimate ##
  if(is.null(phi.prev)){
    beta.old <-  t(x) %*% y / apply(x^2, 2, sum)
    delta.old <- t(xx) %*% y / apply(xx^2, 2, sum)
    bb <- beta.old[jk1] * beta.old[jk2]
    bb[ abs(bb) <= eps1 ] <- eps1
    gamma.old <- delta.old / bb
  }else{
    beta.old <- phi.prev[1:p]
    gamma.old <- phi.prev[(p+1):(p+K)]
    delta.old <- gamma.old * beta.old[jk1] * beta.old[jk2]
  }

  ## Adaptive Weights ##
  if(is.null(phi.weight)){
    phi.weight <- rep(1, p+K)
  }

  beta.w <- phi.weight[1:p]
  delta.w <- phi.weight[(p+1):(p+K)]

  beta.w[ beta.w <= eps1 ] <- eps1

  bb.w <- beta.w[jk1] * beta.w[jk2]

  gamma.w <- delta.w / bb.w
  gamma.w[ gamma.w <= eps1 ] <- eps1


  ## Index sets needed for subfunction entire.beta ##
  NEWIDX1 <- matrix(NA, p-1, p)
  NEWIDX2 <- matrix(NA, K-p+1, p)
  IDX1 <- IDX2 <- matrix(NA, K-p+1, p)
  idx <- seq(1:p)
  for(J in 1:p){
    if(J==1){
        idx1 <- rep(1, p-1)
        idx2 <- c(2:p)
    } else if(J==p){
        idx1 <- c(1:(p-1))
        idx2 <- rep(p, p-1)
    } else{
        idx1 <- c(1:J, rep(J, p-1-J))
        idx2 <- c(rep(J, J-1),c((J+1):p))
    }
    NEWIDX1[,J] <- (idx1-1)*(p-idx1/2) + (idx2-idx1)

    idxJ <- idx[-J]
    i <- 1
    for(j in 1:(p-2)){
      for(k in (j+1):(p-1)){
        IDX1[i,J] <- idxJ[j]
        IDX2[i,J] <- idxJ[k]
        i <- i + 1
      }
    }
    idx1 <- IDX1[,J]
    idx2 <- IDX2[,J]
    NEWIDX2[,J] <- (idx1-1)*(p-idx1/2) + (idx2-idx1)
  }


  ##### Iterative Approximation ######
  iter <- 0
  differ <- 1
  crit.old <- crossprod(y - x%*%beta.old - xx%*%delta.old) + lambda1*sum(abs(beta.old)/beta.w) + lambda2*sum(abs(gamma.old)/gamma.w)

 
  while(differ > eps2){

    ################################
    ### STEP 1: solve for \gamma ###
    ################################
    bb <- beta.old[jk1] * beta.old[jk2]
    bb[ abs(bb) < eps1 ] <- eps1

    ## Construct a new X ##
    newx <- rep(bb, each=n) * xx

    ## Construct a new Y ##
    newy <- y - x %*% beta.old

    ## Calculate a new gamma ##
    gamobj <- solve.gam(newx, newy, lambda2, gamma.w)
    gamma.new <- gamobj$phi.new
    
    ###############################
    ### STEP 2: solve for \beta ###
    ###############################
    beta.new <- beta.old
    beta.new <- sapply(1:p, entire.beta, beta.new, gamma.new, x, xx, y, lambda1, beta.w, NEWIDX1, NEWIDX2, IDX1, IDX2, n)

    delta.new <- gamma.new*beta.new[jk1]*beta.new[jk2]
    crit.new <- crossprod(y - x%*%beta.new - xx%*%delta.new) + lambda1*sum(abs(beta.new)/beta.w) + lambda2*sum(abs(gamma.new)/gamma.w)
    differ <- abs(crit.new - crit.old) / abs(crit.old)

    crit.old <- crit.new
    beta.old <- beta.new
    gamma.old <- gamma.new


    iter <- iter + 1

    if(iter >= max.iter){
      cat("Warning: the whole iteration did not converge.. \n")
      break
    }

  } ## end while

  phi.prev <- c(beta.new, gamma.new)

  gamma.new[abs(gamma.new) <= eps1] <- 0
  beta.new[abs(beta.new) <= eps1] <- 0

  ## Rescale the coeff estimates back into the original scale ##
  beta <- beta.new
  delta <- gamma.new * beta.new[jk1] * beta.new[jk2]

  beta <- beta / normx
  delta <- delta / normxx
  beta0 <- meany - sum(beta*meanx) - sum(delta*meanxx)

  phi <- c(beta0, beta, delta)

  return( list(phi=phi, phi.prev=phi.prev, iter=iter) )
}






##
## subfunctions
##

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


entire.beta <- function(J, beta.new, gamma.new, x, xx, y, lambda1, beta.w, NEWIDX1, NEWIDX2, IDX1, IDX2, n){
      ## calcualte newx ##
      newidx <- NEWIDX1[,J]
      newx <- rep((gamma.new[newidx]*beta.new[-J]), each=n)*xx[,newidx]
      newx <- x[,J] + rowSums(newx)

      ## calculate newy ##
      idx1 <- IDX1[,J]
      idx2 <- IDX2[,J]
      newidx2 <- NEWIDX2[,J]
      bb.mJ <- beta.new[idx1]*beta.new[idx2]
      gamma.mJ <- gamma.new[newidx2]
      xx.mJ <- xx[,newidx2]
      newy <- y - rowSums( rep(beta.new[-J], each=n) * x[,-J] )
      newy <- newy - rowSums( rep((gamma.mJ*bb.mJ), each=n) * xx.mJ )

      beta.newJ <- solve.betaJ(newx, newy, lambda1, beta.w[J])
      return(beta.newJ)
}


solve.betaJ <- function(x, y, lambda, bet){

  xty <- sum(x*y)
  xtx <- sum(x^2)

  if(xty - lambda/(2*bet) > 0){
    phi <- ( xty - lambda/(2*bet) ) / xtx
  }else if(xty + lambda/(2*bet) < 0){
    phi <- ( xty + lambda/(2*bet) ) / xtx
  }else{
    phi <- 0
  }
  return(phi)
}






