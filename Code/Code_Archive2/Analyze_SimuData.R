source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/code/ToolBox_PERCH.R")
library(MVB)
library(doMC)
library(foreach)

registerDoMC(detectCores())

# Simulate Multivariate Bernoulli Distribution
#---------------------------------------------------------------------##
# Dimension of Y, and sample size
K = 7; n = 100
# Dimension of theta
D = 2^K-1
# True value of theta
set.seed(123)
tmp = round(rnorm(D, -(seq(-.5, 1, length.out=D)^2)-1, seq(1,3,length.out=D)),1)
tmp[sample(1:D,floor(D/2))] = 0
theta0 = matrix(tmp,nrow=1)

## design matrix: only has intercept for each parameter
x = matrix(rep(1,1*n),nrow=n)
Lmat = GenMatrices(K)

set.seed(123)
res = mvb.simu(theta0, x, K = K, offset = 0)
dat = res$response
measure = LtoY(dat, TPR=0.9, FPR=0.3)

# True parameter values
phis = exp(Lmat$Lmatrix%*%t(theta0))
TruePi0 = 1/sum(phis)
TruePi1toK = Lmat$PiMat%*%phis*TruePi0
TruePi = round(c(TruePi0, TruePi1toK),3)
TrueMu = as.vector(round(Lmat$MuMat%*%phis*TruePi0,3))
TrueMu
TruePi

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



tid = as.numeric(Sys.getenv("SGE_TASK_ID"))
ParGrid = data.frame(iter=c(rep(1000,4),rep(4500,4)), 
					burnin=c(rep(300,4),rep(500,4)), 
					pi0sigma=rep(c(0.01,0.01,0.07,0.07),2))

# Block-MH sampling with prior on Pi
PiInit = c(0.005,pihat2[1]-0.005, pihat2[2:5], pihat2[6]-0.005, 0.005)
time.start = proc.time()
if (tid %% 2 == 0)
{
	Posterior = post.mu.pi.ByBlock(K=K, mu.init=NULL, pi.init=PiInit , ParMatrix=Lmat, prior.alpha=1.15,
                         iter=ParGrid[tid,1], inner.iter=10, burnin=ParGrid[tid,2], inner.burnin=5, dat=measure, MH.sigmaOfpi0=ParGrid[tid,3])
	time.used = proc.time() - time.start
	print(time.used)	
	print(system("lscpu"))
	filename = paste0("Iter",ParGrid[tid,1],"_measure_",tid,".RData")
} else {
	Posterior = post.mu.pi.ByBlock(K=K, mu.init=NULL, pi.init=PiInit , ParMatrix=Lmat, prior.alpha=1.15,
                         iter=ParGrid[tid,1], inner.iter=10, burnin=ParGrid[tid,2], inner.burnin=5, dat=dat, MH.sigmaOfpi0=ParGrid[tid,3])
	time.used = proc.time() - time.start
	print(time.used)
	print(system("lscpu"))
	filename = paste0("Iter",ParGrid[tid,1],"_dat_",tid,".RData")
}

save(Posterior, measure, dat, time.used, file=filename)



