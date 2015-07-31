source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")
library(MVB)
library(limSolve)
library(doMC)
library(foreach)

# Simulate Multivariate Bernoulli Distribution
#---------------------------------------------------------------------##
# Dimension of Y, and sample size
K = 5; n = 100
# Dimension of theta
D = 2^K-1
# True value of theta
set.seed(123)
Lmat = GenMatrices(K)
unweighted = c(0.1, cumprod(1 - rep(0.5,K)))
pi.prior = round(unweighted/(sum(unweighted)),3)

# compare sampling method
t=proc.time()
Phis1= xsample( E = Lmat$PiMat, F = pi.prior[-1]/pi.prior[1], G = diag(rep(1,Lmat$J1)), H = rep(0,Lmat$J1),
                iter = 500, burnin =50, type = "mirror", test=FALSE)$X
proc.time()-t

t=proc.time()
Phis2 = xsample( E = Lmat$PiMat, F = pi.prior[-1]/pi.prior[1], G = diag(rep(1,Lmat$J1)), H = rep(0,Lmat$J1),
                iter = 2000, burnin = 3000, type = "rda", test=FALSE)$X
proc.time()-t

layout(matrix(1:2,ncol=1))
par(mar=c(2,4,1,1))
hist(Phis1[,2],breaks=15)  

hist(Phis2[,2], breaks=15)
layout(matrix(1,nrow=1))


# calculate canonical parameters
phis = Phis1[2,]
theta = solve(Lmat$Lmatrix, log(phis))
theta0 = theta
#theta0[20:D] = 0

# True parameter values
TruePi0 = 1/(sum(phis)+1)
TruePi1toK = Lmat$PiMat%*%phis*TruePi0
TruePi = round(c(TruePi0, TruePi1toK),3)
TrueMu = as.vector(round(Lmat$MuMat%*%phis*TruePi0,3))
TrueMu
TruePi

## design matrix: only has intercept for each parameter
x = matrix(rep(1,1*n),nrow=n)
res = mvb.simu(matrix(theta0,nrow=1), x, K = K, offset = 0)
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





