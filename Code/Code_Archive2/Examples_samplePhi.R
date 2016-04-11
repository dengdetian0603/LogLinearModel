mu = c(0.8, 0.2, 0.1)
pi0 = 0.15

X0 = rep(0,7)
X0[4] = (1-pi0-mu[3])/pi0
X0[5] = (1-pi0-mu[2])/pi0
X0[6] = (1-pi0-mu[1])/pi0
X0[7] = (sum(mu)+2*pi0-2)/pi0
X0 = matrix(X0)

b1 = matrix(c(1,0,0,-1,-1,0,1))
b2 = matrix(c(0,1,0,-1,0,-1,1))
b3 = matrix(c(0,0,1,0,-1,-1,1))
basis = cbind(b1,b2,b3)

n = 60000

V = runif(n, max=(1-mu[1]-pi0)/pi0)
U = runif(n, max=(1-mu[3]-pi0)/pi0-V)
W.min = (2-2*pi0-sum(mu))/pi0-U-V
W = runif(n, min=W.min, max=(1-mu[2]-pi0)/pi0-V-U)

UVW = rbind(U,V,W)[,!is.na(W)&!is.na(U)]
N = ncol(UVW)

candidate = apply(basis%*%UVW,2,function(x) x + X0)
index = apply(candidate,2,function(x) sum(x>0&x<(1/pi0-1))>6)
phi = candidate[,index]

layout(matrix(1:9,nrow=3,byrow=TRUE))
hist(phi[1,],breaks=50)
hist(phi[2,],breaks=50)
hist(phi[3,],breaks=50)
hist(phi[4,],breaks=50)
hist(phi[5,],breaks=50)
hist(phi[6,],breaks=50)
hist(phi[7,],breaks=50)
layout(matrix(1,nrow=1))

## plot the feasible space of (u,v,w)
library(rgl)
open3d()
plot3d(x=UVW[1,index],y=UVW[2,index],z=UVW[3,index],
       xlab = "u", ylab = "v", zlab = "w", 
       col="steelblue", type="s", size=1,alpha=0.75)

pairs(t(candidate[,index]),pch=20)

## plot thetas generated from uniform (u,v,w)
#library(pracma)
Lmat = matrix(c(1,0,0,0,0,0,0,
                0,1,0,0,0,0,0,
                0,0,1,0,0,0,0,
                1,1,0,1,0,0,0,
                1,0,1,0,1,0,0,
                0,1,1,0,0,1,0,
                1,1,1,1,1,1,1), nrow=7,byrow=TRUE)

thetas = apply(phi,2,function(x) solve(Lmat,log(x)))
plot3d(x=thetas[1,],y=thetas[2,],z=thetas[3,],
       xlab = "theta1", ylab = "theta2", zlab = "theta3", 
       col="steelblue", type="s", size=1,alpha=0.75)

# ----------------------------------------------------------------------------------- #

library(limSolve)
library(doMC)
library(foreach)
source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")

registerDoMC(4) 

K = 3 # 
mu.prior = c(0.79,0.79,0.79,0.79,0.79)
pi.prior = c(0.2,0.5,0.299,0.001)


#K = 2 # when K < 3, only single solution
#mu.prior = c(0.4,0.3)
#pi.prior = c(0.4,0.5,0.1)

#-- get all possible combinations
L.all = as.matrix(AllComb(K))[-1,]
S = apply(L.all,1,sum)
index = order(S)
L.all = L.all[index,]
S = S[index]
J1 = length(S)
      
#-- expand all possible combinations
L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
      expandL.LL(L.all[i,])
}

PiMat = matrix(NA,nrow=K,ncol=J1)
for (i in 1:K)
{
      PiMat[i,] = as.numeric(S==i)
}

MuMat = matrix(NA,nrow=K,ncol=J1)
for (i in 1:K)
{
      MuMat[i,] = as.numeric(L.all.exp[,i]==1)
}

## Given Pi only
Phis = xsample( E = PiMat, F = pi.prior[-1]/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                iter = 80, output = 40, type = "cda")
mu.sample = MuMat%*%t(Phis$X)*pi.prior[1]
boxplot(t(mu.sample))

## Given Mu and Pi0 only
# need additional constrain that pis sum up to 1
mu.prior = mu.sample[,1]
Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi.prior[1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                iter = 60, output = 40, type = "cda")
pi.sample = PiMat%*%t(Phis$X)*pi.prior[1]
round(pi.sample,3)
boxplot(t(pi.sample))

## Given both Mu and Pi
mu.prior = mu.sample[,2]
Phis = tryCatch(xsample( E = rbind(PiMat,MuMat), F = c(pi.prior[-1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
               iter = 300, output = 10, type = "cda"),error=function(e) "incompatible constraints")
Phis$X

tryCatch(xsample( E = rbind(PiMat,MuMat), F = c(0.2148334, 0.2990931, 0.1139873, 0.5923605, 0.4486338, 0.1139873)/0.3720862, G = diag(rep(1,J1)), H = rep(0,J1),
         iter = 300, output = 10, type = "cda") ,error=function(e) "incompatible constraints")



L.all

i=10
cell.prob = c(1, Phis$X[i,])*pi.prior[1]
rmultinom(1,100,cell.prob)

apply(L.all,1,function(x) paste(as.character(x),collapse=""))

rDataWithPiMu.LL(5, mu.prior, pi.prior,1)

## solve for qe prameters
library(nloptr)

# Objective function
eval_f = function(x, xx, B, Lmat, b) 
{
      list("objective" = sum((x-xx)^2),
           "gradient" = 2*(x-xx))
}

# non-linear equality constraints
eval_g_eq = function(x, B, Lmat, b,xx)
{
      constr = B%*%exp(Lmat%*%x) - b
      grad = B%*%diag(exp(Lmat%*%x)[,1])%*%Lmat
      return( list( "constraints"=constr, "jacobian"=grad ) )
}

# initial values
npar = (K*(K+1)/2)
x0 = Phis$X[2,1:npar]
# lower and upper bounds of control
lb = rep(-Inf,npar)
ub = rep(Inf,npar)
local_opts = list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-5 )
opts = list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-5,
              "maxeval" = 20000,
              "local_opts" = local_opts )
res = nloptr( x0=x0, eval_f=eval_f,lb=lb, ub=ub,
              eval_g_eq=eval_g_eq, opts=opts,
              xx=x0, B=rbind(PiMat,MuMat), b=c(pi.prior[-1],mu.prior)/pi.prior[1], Lmat = L.all.exp[,1:npar])
res


# -------------------------------------------------------------------------------- #
K = 3 # 
mu.prior = mu
pi.prior = pi0

#-- get all possible combinations
L.all = as.matrix(AllComb(K))[-1,]
S = apply(L.all,1,sum)
index = order(S)
L.all = L.all[index,]
S = S[index]
J1 = length(S)

#-- expand all possible combinations
L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
      expandL.LL(L.all[i,])
}

PiMat = matrix(NA,nrow=K,ncol=J1)
for (i in 1:K)
{
      PiMat[i,] = as.numeric(S==i)
}

MuMat = matrix(NA,nrow=K,ncol=J1)
for (i in 1:K)
{
      MuMat[i,] = as.numeric(L.all.exp[,i]==1)
}


## Given Mu and pi0 only
# need additional constrain that pis sum up to 1
Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi.prior[1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                iter = 6000, output = 1000, type = "mirror")
pi.sample = PiMat%*%t(Phis$X[1:1000,])*pi.prior[1]
boxplot(t(pi.sample))
mu.sample = MuMat%*%t(Phis$X[1:1000,])*pi.prior[1]
boxplot(t(mu.sample))

## plot thetas generated from uniform phi
thetas = apply(Phis$X,1,function(x) solve(L.all.exp,log(x)))
plot3d(x=thetas[1,],y=thetas[2,],z=thetas[3,],
       xlab = "theta1", ylab = "theta2", zlab = "theta3", 
       col="steelblue", type="s", size=1,alpha=0.75)

plot3d(x=thetas[4,],y=thetas[5,],z=thetas[6,],
       xlab = "theta4", ylab = "theta5", zlab = "theta6", 
       col="steelblue", type="s", size=1,alpha=0.75)

layout(matrix(1:9,nrow=3,byrow=TRUE))
hist(thetas[1,],breaks=50)
hist(thetas[2,],breaks=50)
hist(thetas[3,],breaks=50)
hist(thetas[4,],breaks=50)
hist(thetas[5,],breaks=50)
hist(thetas[6,],breaks=50)
hist(thetas[7,],breaks=50)
layout(matrix(1,nrow=1))

GetVol.mupi0(c(0.8,0.2,0.1),0.15,2000)