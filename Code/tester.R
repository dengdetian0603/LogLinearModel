library(dplyr)
library(magrittr)
library(gtools)
library(doMC)
library(foreach)
library(MVB)
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/LS.model.R")
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/LL.exch_150217.R")
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/StateToMeasure.R")
registerDoMC(4) 

Pis = c(0.1,0.5,0.3,0.05,0.04,0.009,0.001-1e-5,1e-5)
K = length(Pis)-1

n = 500; m = 50000

# sim strategy 1
L = mvb.exch.sim(Pis, n)
# sim strategy 2
theta0 = matrix(c(rnorm(K,-0.3,1),rnorm(choose(K,2),-0.5,1.5)),nrow=1)
x1 = matrix(rep(1,1*n),nrow=n)
L = mvb.simu(theta0, x1, K, offset = 0)$response

table(factor(apply(L,1,sum),levels=0:K))

tpr = pmin(rnorm(K,0.9,0.1),1)
fpr = pmax(rnorm(K,0.1,0.05),0)
Y = LtoY(L,tpr,fpr) # generating Y with TPR and FPR 
L=Y

#-- expand the observed vectors
L.expand = foreach(i=1:nrow(L), .combine=rbind) %dopar% {
      expandL(L[i,])
}
#-- get all possible combinations
L.all = as.matrix(AllComb(ncol(L))) 
#-- expand all possible combinations
L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
      expandL(L.all[i,])
}

#--------------------- Log likelihood with L -------------------------#
fn_L = function(par){loglik1(par,L=L.expand, L.All=L.all.exp)}

fit1 = optim(par=rep(0,K+choose(K,2)), fn=fn_L, method="BFGS",control=list(maxit=500,abstol=1e-10))
beta1 = fit1$par

#fitMVB = mvbfit(x1, L, output = 1);fitMVB$beta

# Expected marginal total S
true.pi = QE.ThetaToPi(theta0,L.all.exp,K)
S.freq = round(m*true.pi,0)

#--------------------- Log likelihood with L and S ---------------------#
fn_LS = function(par) {loglik2(par, L=L.expand, L.All=L.all.exp, S=S.freq)}
fit2 = optim(par=rep(0,K+choose(K,2)), fn=fn_LS, method="BFGS",control=list(maxit=500,abstol=1e-10))
beta2 = fit2$par



plot(beta1,ylim=c(min(beta1,beta2,theta0[1,]),max(beta1,beta2,theta0[1,]))) # black: L model
points(beta2,col="blue") # blue: L + S model
points(as.vector(theta0),col="red") # red: true value

sum((beta1-theta0[1,])^2);sum((beta2-theta0[1,])^2)
