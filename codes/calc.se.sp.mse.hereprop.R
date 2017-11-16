setwd('D:/KC_file/Dropbox/VSSH/results/cond8')

# for proportion of violation #
source(paste(dirname(dirname(getwd())),"/codes/calc.prop.heredity.maineffect.r",sep=""))
# among num.sim numbers of simulation #
# how many times heredity is violated #

set.seed(1)

q = 10
numsig = 4  # number of significant main effects
numnonusemain = 0  # number of mains that dont have interactions

main.magnitude = 4
int.magnitude = 1
quad.magnitude = 1

beta1 = c(rep(1,numsig),rep(0,q-numsig))
beta2 = combn(c(rep(1,numsig-numnonusemain),rep(0,q-numsig+numnonusemain)),2,prod)
beta3 = c(rep(1,numsig),rep(0,q-numsig))
beta1 = beta1 * main.magnitude
beta2 = beta2 * int.magnitude
beta3 = beta3 * quad.magnitude
beta = c(beta1, beta2, beta3)

sim.num = 50

coef.mat.shim.poly = t(read.csv("coef.plain.shim.csv",header=T)[,-1])
coef.mat.lassoori.poly = t(read.csv("coef.plain.lasso.csv",header=T)[,-1])
coef.mat.lasso.poly = t(read.csv("coef.vanilla.lasso.csv",header=T)[,-1])
coef.mat.step.poly = t(read.csv("coef.vanilla.stepwise.csv",header=T)[,-1])
coef.mat.stepori = t(read.csv("coef.plain.stepwise.csv",header=T)[,-1])

mse = read.csv("mse.csv",header=T)[,-1]

eps = c(0,0,0)
#eps = c(main.magnitude, int.magnitude, quad.magnitude)/100

se = matrix(,5,sim.num)
sp = matrix(,5,sim.num)
heredity.violation = matrix(,5,sim.num)
# each column is one iteration
# each row is one method
coef = list(coef.mat.lasso.poly,coef.mat.lassoori.poly,coef.mat.step.poly,coef.mat.stepori,coef.mat.shim.poly)
for (i in 1:5) {
se[i,] = (apply(abs(coef[[i]][,1:q][,which(beta1!=0)]) > eps[1],1,sum)+
	apply(abs(coef[[i]][,(q+1):(q*(q+1)/2)][,which(beta2!=0)]) > eps[2],1,sum)+
	apply(abs(coef[[i]][,(q*(q+1)/2+1):(q*(q+3)/2)][,which(beta3!=0)]) > eps[3],1,sum))/sum(beta!=0)
sp[i,] = (apply(abs(coef[[i]][,1:q][,which(beta1==0)]) <= eps[1],1,sum)+
	apply(abs(coef[[i]][,(q+1):(q*(q+1)/2)][,which(beta2==0)]) <= eps[2],1,sum)+
	apply(abs(coef[[i]][,(q*(q+1)/2+1):(q*(q+3)/2)][,which(beta3==0)]) <= eps[3],1,sum))/sum(beta==0)
heredity.violation[i,] = proportion_heredity(coef[[i]],q)
}
rownames(se) = c("se.vanilla_lasso", "se.plain_lasso", "se.vanilla_stepwise", "se.plain_stepwise", "se.SHIM")
rownames(sp) = c("sp.vanilla_lasso", "sp.plain_lasso", "sp.vanilla_stepwise", "sp.plain_stepwise", "sp.SHIM")
rownames(heredity.violation) = c("vanilla_lasso", "plain_lasso", "vanilla_stepwise", "plain_stepwise", "SHIM")

se.table = rbind(apply(se,1,mean),(apply(se,1,sd)/sqrt(sim.num)))
rownames(se.table) = c("mean", "s.e.")
sp.table = rbind(apply(sp,1,mean),(apply(sp,1,sd)/sqrt(sim.num)))
rownames(sp.table) = c("mean", "s.e.")
write.csv(round(t(cbind(se.table,sp.table)),3),"sesp.csv")


mse.table = rbind(apply(mse,2,mean),(apply(mse,2,sd)/sqrt(50)))
rownames(mse.table) = c("mean", "s.e.")
colnames(mse.table) = c("vanilla_lasso", "plain_lasso", "vanilla_stepwise", "plain_stepwise", "SHIM")
write.csv(round(t(mse.table),2),"mse.finaltable.csv")

write.csv(1-t(heredity.violation),"heredity.violation.csv")
here.table = apply(heredity.violation,1,mean)
write.csv(1-here.table,"here.maintain.final.csv")

## proportion_heredity(coef[[2]],q)


################### NOT USE ####################
#se.shim.poly.m = apply(abs(coef.mat.shim.poly[,which(beta[1:q]!=0)])>cnt.criterion,1,sum) /sum(beta[1:q]!=0)
#sp.shim.poly.m = apply(abs(coef.mat.shim.poly[,which(beta[1:q]==0)])<=cnt.criterion,1,sum) /sum(beta[1:q]==0)
#se.ourlasso.poly.m = apply(abs(coef.mat.lasso.poly[,which(beta[1:q]!=0)])>cnt.criterion,1,sum) /sum(beta[1:q]!=0)
#sp.ourlasso.poly.m = apply(abs(coef.mat.lasso.poly[,which(beta[1:q]==0)])<=cnt.criterion,1,sum) /sum(beta[1:q]==0)
#se.ourstep.poly.m = apply(abs(coef.mat.step.poly[,which(beta[1:q]!=0)])>cnt.criterion,1,sum) /sum(beta[1:q]!=0)
#sp.ourstep.poly.m = apply(abs(coef.mat.step.poly[,which(beta[1:q]==0)])<=cnt.criterion,1,sum) /sum(beta[1:q]==0)
#se.lasso.poly.m = apply(abs(coef.mat.lassoori.poly[,which(beta[1:q]!=0)])>cnt.criterion,1,sum) /sum(beta[1:q]!=0)
#sp.lasso.poly.m = apply(abs(coef.mat.lassoori.poly[,which(beta[1:q]==0)])<=cnt.criterion,1,sum) /sum(beta[1:q]==0)
#se.step.poly.m = apply(abs(coef.mat.stepori[,which(beta[1:q]!=0)])>cnt.criterion,1,sum) /sum(beta[1:q]!=0)
#sp.step.poly.m = apply(abs(coef.mat.stepori[,which(beta[1:q]==0)])<=cnt.criterion,1,sum) /sum(beta[1:q]==0)

#se.shim.poly.2fi = apply(abs(coef.mat.shim.poly[,which(beta[(q+1):(q*(q+1)/2)]!=0)+q])>cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]!=0)
#sp.shim.poly.2fi = apply(abs(coef.mat.shim.poly[,which(beta[(q+1):(q*(q+1)/2)]==0)+q])<=cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]==0)
#se.ourlasso.poly.2fi = apply(abs(coef.mat.lasso.poly[,which(beta[(q+1):(q*(q+1)/2)]!=0)+q])>cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]!=0)
#sp.ourlasso.poly.2fi = apply(abs(coef.mat.lasso.poly[,which(beta[(q+1):(q*(q+1)/2)]==0)+q])<=cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]==0)
#se.ourstep.poly.2fi = apply(abs(coef.mat.step.poly[,which(beta[(q+1):(q*(q+1)/2)]!=0)+q])>cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]!=0)
#sp.ourstep.poly.2fi = apply(abs(coef.mat.step.poly[,which(beta[(q+1):(q*(q+1)/2)]==0)+q])<=cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]==0)
#se.lasso.poly.2fi = apply(abs(coef.mat.lassoori.poly[,which(beta[(q+1):(q*(q+1)/2)]!=0)+q])>cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]!=0)
#sp.lasso.poly.2fi = apply(abs(coef.mat.lassoori.poly[,which(beta[(q+1):(q*(q+1)/2)]==0)+q])<=cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]==0)
#se.step.poly.2fi = apply(abs(coef.mat.stepori[,which(beta[(q+1):(q*(q+1)/2)]!=0)+q])>cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]!=0)
#sp.step.poly.2fi = apply(abs(coef.mat.stepori[,which(beta[(q+1):(q*(q+1)/2)]==0)+q])<=cnt.criterion,1,sum) /sum(beta[(q+1):(q*(q+1)/2)]==0)

#se.shim.poly.q = apply(abs(coef.mat.shim.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)+q*(q+1)/2])>cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)
#sp.shim.poly.q = apply(abs(coef.mat.shim.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)+q*(q+1)/2])<=cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)
#se.ourlasso.poly.q = apply(abs(coef.mat.lasso.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)+q*(q+1)/2])>cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)
#sp.ourlasso.poly.q = apply(abs(coef.mat.lasso.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)+q*(q+1)/2])<=cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)
#se.ourstep.poly.q = apply(abs(coef.mat.step.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)+q*(q+1)/2])>cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)
#sp.ourstep.poly.q = apply(abs(coef.mat.step.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)+q*(q+1)/2])<=cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)
#se.lasso.poly.q = apply(abs(coef.mat.lassoori.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)+q*(q+1)/2])>cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)
#sp.lasso.poly.q = apply(abs(coef.mat.lassoori.poly[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)+q*(q+1)/2])<=cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)
#se.step.poly.q = apply(abs(coef.mat.stepori[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)+q*(q+1)/2])>cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]!=0)
#sp.step.poly.q = apply(abs(coef.mat.stepori[,which(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)+q*(q+1)/2])<=cnt.criterion,1,sum) /sum(beta[(q*(q+1)/2+1):(q*(q+3)/2)]==0)
################### NOT USE ####################


