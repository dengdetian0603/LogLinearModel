source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")
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

# MLE of saturated model on perfect data
muhat = apply(dat, 2, mean)
pihat = table(apply(dat,1,sum))/n
pihat[(length(pihat)+1):(K+1)] = 0

# MLE on imperfect data
muhat2 = apply(measure, 2, mean)
pihat2 = table(apply(measure,1,sum))/n
pihat2[(length(pihat2)+1):(K+1)] = 0

# prior
Sigma = 1.6
Alpha = 1.14
x = rnorm(1000*K, 0, Sigma)
mu = matrix(exp(x)/(1+exp(x)),ncol=K)
pi0 = apply(mu, 1, function(mu) runif(1,0,1-max(mu)))
pi1toK = sapply(pi0, function(x) {tmp=rStickBreak(K-1, Alpha); 
                                  tmp2=c(tmp,1-sum(tmp)); 
                                  return(tmp2*(1-x))})
prior = apply(cbind(mu,pi0,t(pi1toK)),2,mean)

#---------------------- Chain Diagonosis -------------------------------- # 
load("~/Downloads/Iter1000_dat_1.RData")
tmp = Posterior
apply(tmp$history.accept[1:700,],2,mean)

library(coda)
mu.chain = mcmc(tmp$posterior[,1:15])
summary(mu.chain)
plot(mu.chain)
autocorr.plot(mu.chain)
rejectionRate(mu.chain)

mc.est = round(summary(mu.chain)$statistics[,1],3)
est.compare = cbind(c(TrueMu, TruePi), c(muhat, pihat), mc.est, round(prior,3), c(muhat2,pihat2))
rownames(est.compare) = names(mc.est)
colnames(est.compare) = c("true","mle.gold","mc.est","prior","mle.brown")
est.compare


# chain1 =  mcmc(tmp$posterior[1:2000,5])
# chain2 =  mcmc(tmp$posterior[-(1:2000),5])
# combinedchains = mcmc.list(chain1, chain2) 
# plot(combinedchains)
# gelman.diag(combinedchains)
# gelman.plot(combinedchains)

