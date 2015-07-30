## 
library(limSolve)
library(doMC)
library(foreach)
source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")

registerDoMC(4) 

dat = matrix(c(rep(c(1,0,0),10),rep(c(0,1,0),5),rep(c(0,0,1),2), rep(0,9),
               rep(c(1,1,0),4), rep(c(1,0,1),2), rep(c(0,1,1),1)), byrow=TRUE, ncol=3)
muhat = apply(dat,2,mean)

mat = GenMatrices(3)

density.YgivenMuPi(K=3, dat=dat, y=NULL, mu=mu.prior, pi=pi.prior, logscale=FALSE, burnin=300, n=100, method="cda", parMat=mat)
density.YgivenMuPi(K=3, dat=dat, y=NULL, mu=muhat, pi=c(3,17,7,0.1)/27.1, logscale=FALSE, burnin=300, n=100, method="cda", parMat=mat)

dStickBreak(pi.prior[-1]/(1-pi.prior[1]), 1.1)
prior.MuPi(mu=mu.prior, pi=pi.prior, Sigma=rep(1.6,3), Alpha=1.1, logscale=FALSE)

density.YMuPi(K=3, dat=dat, mu=mu.prior,pi=pi.prior, SigmaInPrior=rep(1.6,3), AlphaInPrior=1.1,
               logscale=TRUE, inner.burnin=200, inner.iter=100, method="cda", ParMat=mat)


# ---------------- Build a simple sampler -------------------------- #
source("~/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")

# TODO:
#tmp = post.mu.pi(K=3, mu.init=NULL, pi.init=c(3,17,7,0.1)/27.1, ParMatrix=mat, prior.alpha=1.15,
#                  iter=45, inner.iter=10, burnin=5, inner.burnin=5, dat=dat, MH.par=c(1.3,1.3,1.6,0.15))
tmp = post.mu.pi.ByBlock(K=3, mu.init=NULL, pi.init=c(3,17,7,0.1)/27.1, ParMatrix=mat, prior.alpha=1.15,
                 iter=5000, inner.iter=10, burnin=1000, inner.burnin=5, dat=dat, MH.par=c(1.3,1.3,1.6,0.14))

# -------------------- Chain visualization and diagonosis ----------------------------------- #
#tmp = postlist[[18]]
load("~/Downloads/Iter1000_dat_1.RData")
tmp = Posterior
apply(tmp$history.accept[1:700,],2,mean)

layout(matrix(1:2,ncol=1))
plot(density(tmp$posterior[1:700,1],adjust=1.8), xlab="value of parameter", xlim=c(0,1),main="Posterior Density", ylim=c(0,7.5))
lines(density(tmp$posterior[1:700,2],adjust=1.8),col=2)
lines(density(tmp$posterior[1:700,3],adjust=1.8),col=4)
legend("topright",legend=c(paste0("Mu1 : ",round(mean(tmp$posterior[1:700,1]),3)),
                           paste0("Mu2 : ",round(mean(tmp$posterior[1:700,2]),3)),
                           paste0("Mu3 : ",round(mean(tmp$posterior[1:700,3]),3))),
       lty=1, col=c(1,2,4))
plot(density(tmp$posterior[1:700,8],adjust=1.8), xlab="value of parameter", xlim=c(0,1), ylim=c(0,15), main="Posterior Density")
lines(density(tmp$posterior[1:700,9],adjust=1.8),col=2)
lines(density(tmp$posterior[1:700,10],adjust=1.8),col=4)
lines(density(tmp$posterior[1:700,11],adjust=1.8),col=3)
legend("topright",legend=c(paste0("Pi0 : ",round(mean(tmp$posterior[1:700,8]),3)),
                           paste0("Pi1 : ",round(mean(tmp$posterior[1:700,9]),3)),
                           paste0("Pi2 : ",round(mean(tmp$posterior[1:700,10]),3)),
                           paste0("Pi3 : ",round(mean(tmp$posterior[1:700,11]),3))),
       lty=1, col=c(1,2,4,3))
layout(matrix(1))


layout(matrix(1:4,ncol=1))
par(mar=c(2,4,1,1))
plot(tmp$posterior[1:700,1],type="l",xlab="iteration",ylab="Mu1")
plot(tmp$posterior[1:700,2],type="l",xlab="iteration",ylab="Mu2")
plot(tmp$posterior[1:700,4],type="l",xlab="iteration",ylab="Pi0")
plot(tmp$posterior[1:700,5],type="l",xlab="iteration",ylab="Pi1")
layout(matrix(1))


library(coda)
mu.chain = mcmc(tmp$posterior[,1:15])
summary(mu.chain)
plot(mu.chain)
autocorr.plot(mu.chain)
rejectionRate(mu.chain)

chain1 =  mcmc(tmp$posterior[1:2000,5])
chain2 =  mcmc(tmp$posterior[-(1:2000),5])
combinedchains = mcmc.list(chain1, chain2) 
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)





