setwd('D:/Schooling/Higher Education/!papers&research/SHIM_New')

## x2[,1]:x1*x2, x2[,2]:x1*x3, x2[,3]:x1*x4
## x2[,10:11]: x2*x3, x2*x4
## x2[,18]: x3*x4

## mycase2
beta1 = c(1,1,1,1,rep(0,6))*10  # may change
beta2 = rep(0,45)
beta2[c(1,2,3,10,11,18)] = rep(10,6)
beta = c(beta1, beta2)

q = 10

info.4 = as.matrix(read.table('infomat.mycase2.txt'))
coef.4 = as.matrix(read.table('coefmat.mycase2.txt'))

info.lasso.4 = as.matrix(read.table('infomat.lasso.mycase2.txt'))
coef.lasso.4 = as.matrix(read.table('coefmat.lasso.mycase2.txt'))

bmat = coef.4[,-1]
apply(abs(bmat)[,abs(beta)>0]>0,2,sum)
apply(abs(bmat)[,abs(beta)==0]>0,2,sum)

bmat = coef.4[,-1]
bmat2 = coef.lasso.4[,-1]

## get true nonzero for our methods
bmat.main = bmat[,1:q]
bmat.int = bmat[,-(1:q)]

bmat1 = bmat
for (j in 1:q) {
	tmp1 = rep(1,q)
	tmp1[j] = 0
	tmp1 = matrix(tmp1,nrow=1)
	tmp2 = t(apply(tmp1, 1, combn, 2, prod))
	bmat1[,j] = (apply(abs(bmat.int[,tmp2==0])>0,1,sum)>0)*T
}
bmat1[,1:q] = bmat1[,1:q] + (abs(bmat[,1:q])>0)*T

s1 = c(mean(apply(abs(bmat1[,abs(beta)>0])>0,1,sum)), sd(apply(abs(bmat1[,abs(beta)>0])>0,1,sum)), mean(apply(abs(bmat1[,abs(beta)==0])==0,1,sum)), sd(apply(abs(bmat1[,abs(beta)==0])==0,1,sum)))
s2 = c(mean(apply(abs(bmat2[,abs(beta)>0])>0,1,sum)), sd(apply(abs(bmat2[,abs(beta)>0])>0,1,sum)), mean(apply(abs(bmat2[,abs(beta)==0])==0,1,sum)), sd(apply(abs(bmat2[,abs(beta)==0])==0,1,sum)))
cbind(s1, s2)

> cbind(s1, s2)
            s1         s2
[1,]  7.380000  7.0500000
[2,]  1.489221  0.9679198
[3,] 31.990000 32.9400000
[4,]  4.461949  4.3410654

