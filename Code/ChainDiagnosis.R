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
load("~/Documents/workspace/Iter4500_dat_5.RData")
tmp = Posterior
apply(tmp$history.accept[3000:4000,],2,mean)

library(coda)
mu.chain = mcmc(tmp$posterior[3000:4000,])
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

