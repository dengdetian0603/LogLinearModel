source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")
library(MVB)
library(limSolve)
library(doMC)
library(foreach)

# Simulate Multivariate Bernoulli Distribution
#---------------------------------------------------------------------##
# Dimension of Y, and sample size
K = 7; n = 100
# Dimension of theta
D = 2^K-1
# True value of theta
set.seed(123)
Lmat = GenMatrices(K)
unweighted = c(0.1, cumprod(1 - rep(0.5,K)))
pi.prior = round(unweighted/(sum(unweighted)),3)
Phis = xsample( E = Lmat$PiMat, F = pi.prior[-1]/pi.prior[1], G = diag(rep(1,Lmat$J1)), H = rep(0,Lmat$J1),
                iter = 100, output = 100, type = "mirror")$X
phis = Phis[73,]
# True parameter values
TruePi0 = 1/(sum(phis)+1)
TruePi1toK = Lmat$PiMat%*%phis*TruePi0
TruePi = round(c(TruePi0, TruePi1toK),3)
TrueMu = as.vector(round(Lmat$MuMat%*%phis*TruePi0,3))
TrueMu
TruePi

## design matrix: only has intercept for each parameter
x = matrix(rep(1,1*n),nrow=n)
res = mvb.simu(theta0, x, K = K, offset = 0)
dat = res$response
measure = LtoY(dat, TPR=0.9, FPR=0.3)

# MLE of saturated model on perfect data
muhat = apply(dat, 2, mean)
pihat = table(apply(dat,1,sum))/n
muhat
pihat

# MLE on imperfect data
muhat2 = apply(measure, 2, mean)
pihat2 = table(apply(measure,1,sum))/n
muhat2
pihat2

# Block-MH sampling with prior on Pi
PiInit = c(0.005,pihat2[1]-0.005, pihat2[2:5], pihat2[6]-0.005, 0.005)
tmp = post.mu.pi.ByBlock(K=K, mu.init=NULL, pi.init=PiInit , ParMatrix=Lmat, prior.alpha=1.15,
                         iter=7000, inner.iter=10, burnin=3000, inner.burnin=5, dat=measure, MH.par=c(rep(1.3,K),0.14))





