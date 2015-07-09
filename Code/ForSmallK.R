mu = c(0.8, 0.2, 0.1)
pi0 = 0.001

X0 = rep(0,7)
X0[1] = (1-pi0-mu[2]-mu[3])/pi0
X0[2] = (1-pi0-mu[1]-mu[3])/pi0
X0[3] = mu[3]/pi0
X0[4] = (sum(mu)-1+pi0)/pi0
X0 = matrix(X0)

b1 = matrix(c(0,1,-1,-1,1,0,0))
b2 = matrix(c(1,0,-1,-1,0,1,0))
b3 = matrix(c(1,1,-1,-2,0,0,1))
basis = cbind(b1,b2,b3)

n = 20000

V = runif(n, max=(1-mu[1]-pi0)/pi0)
U = runif(n, max=mu[3]/pi0-V)
W.min = apply(rbind((mu[2]+mu[3]+pi0-1)/pi0-V,(mu[1]+mu[3]+pi0-1)/pi0-U),2,max)
W = runif(n, min=W.min, max=mu[3]/pi0-V-U)

UVW = rbind(U,V,W)[,!is.na(W)&!is.na(U)]

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


## 
library(limSolve)
library(doMC)
library(foreach)
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")

registerDoMC(4) 

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



#####

Bmat = matrix(c(1,1,1,1,1,1,1,
                1,0,1,1,0,0,1,
                1,1,0,1,0,1,0,
                1,1,1,0,1,0,0), nrow=4, byrow=TRUE)
rref(Bmat)


