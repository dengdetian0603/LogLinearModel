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

GetVol.mupi0(c(0.8,0.2,0.1),0.15,2000)

dat = matrix(c(rep(c(1,0,0),10),rep(c(0,1,0),5),rep(c(0,0,1),2), rep(0,9),
               rep(c(1,1,0),4), rep(c(1,0,1),2), rep(c(0,1,1),1)), byrow=TRUE, ncol=3)
apply(dat,2,mean)

dYgivenMuPi0(K=3,dat=dat,mu=c(0.59,0.37,0.19),pi0=0.11,n=2000,logscale=TRUE)

# ---------------- Build a simple sampler -------------------------- #
density.YMuPi0(K=3, dat=dat, mu=c(0.49,0.27,0.59),pi0=0.09, Sigma=rep(1.6,3),
               logscale=TRUE, inner.burnin=1000, inner.iter=500, method="cda")

source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")


tmp = post.mu.pi0(K=3, mu.init=c(0.49,0.27,0.59), pi0.init=0.15, 
            iter=1500, inner.iter=500, burnin=500, inner.burnin=100, dat=dat, MH.par=c(1.3,1.3,1.6,0.15))

tmp = postlist[[18]]
apply(tmp$history.accept,2,mean)

plot(density(tmp$posterior[,1],adjust=1.5), xlab="value of parameter", xlim=c(0,1),main="Posterior Density", ylim=c(0,7.5))
lines(density(tmp$posterior[,2],adjust=1.5),col=2)
lines(density(tmp$posterior[,3],adjust=1.5),col=4)
lines(density(tmp$posterior[,4],adjust=1.5),col=3)
legend("topright",legend=c(paste0("Mu1 : ",round(mean(tmp$posterior[,1]),3)),
                           paste0("Mu2 : ",round(mean(tmp$posterior[,2]),3)),
                           paste0("Mu3 : ",round(mean(tmp$posterior[,3]),3)),
                           paste0("Pi0 : ",round(mean(tmp$posterior[,4]),4))),
       lty=1, col=c(1,2,4,3))

layout(matrix(1:4,ncol=1))
par(mar=c(2,4,1,1))
plot(tmp$posterior[,1],type="l",xlab="iteration",ylab="Mu1")
plot(tmp$posterior[,2],type="l",xlab="iteration",ylab="Mu2")
plot(tmp$posterior[,3],type="l",xlab="iteration",ylab="Mu3")
plot(tmp$posterior[,4],type="l",xlab="iteration",ylab="Pi0")
layout(matrix(1))
