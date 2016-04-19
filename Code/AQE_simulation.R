library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
library(nleqslv)
library(BB)
library(mgcv)

sourceCpp("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_Components.cpp")
source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_models_conditional.R")
#---------------------------------------------------------------------------
K=5; Smax=4
dmat = AQE.DesignMatrix(5,4)
cellLabel = getCellLabel(K,dmat$Lmat)

ss_tpr = c(0.05, 0.12, 0.08, 0.15, 0.1); bs_tpr = c(0.8,0.6,0.7,0.7,0.5); bs_fpr = c(0.5, 0.55, 0.4, 0.35, 0.45)
#theta2 = rep(0,10)
theta2=c(-0.1,-0.3, -0.4, -0.8, -0.45, -0.3, -0.4, -0.7, 0.4, -0.9)

Data = AQE.simulate(K=5, Smax=5, N=1200, theta2=theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=134)

PAR0 = list(beta=as.vector(Data$Beta), theta2 = theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr)
HYPAR = list(a=rep(1,K), b=rep(10,K), c=rep(2,K), d=rep(1,K), e=rep(1,K), f=rep(1,K), 
             varbeta=20, vartheta=15, mutheta=-0.1, pmix=0.7)

ctrlpar = data.frame(burnin=1, iter=1000, epsilon=0.002, Nstep=3)
tmp = GetPosterior(K=5, Smax=4, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)

ctrlpar = data.frame(maxiter=100, tol=1e-5)
tmp = GetEMsolution(K=5, Smax=4, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)

# ----------------------------------------------------------------------------

LUmat = cbind(dmat$Lmat, dmat$Umat)
MuMat = dmat$MuMat
Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
cellLabel = apply(Lmat.withZero,1,function(x) paste(as.character(x),collapse=""))

MSS = Data$MSS.case
MBS = Data$MBS.case
MBS.ctrl = Data$MBS.ctrl
X = Data$X
X.unique = uniquecombs(X)
X.index = attr(X.unique, "index")

D = ncol(X)
N.case = nrow(MBS)
N.ctrl = nrow(MBS.ctrl)

maxiter = 100
tol = 1e-5

HyperPar = HYPAR

aa = HyperPar$a; bb = HyperPar$b; cc = HyperPar$c;
dd = HyperPar$d; ee= HyperPar$e; ff = HyperPar$f;
varbeta=HyperPar$varbeta; vartheta=HyperPar$vartheta 
mutheta=HyperPar$mutheta; pmix=HyperPar$pmix


ParInit = PAR0
THETA = c(ParInit$beta, ParInit$theta2, ParInit$ss_tpr, ParInit$bs_tpr, ParInit$bs_fpr)
Beta_index = 1:(K*D)
theta2_index = (1:choose(K,2)) + (K*D)
ss_tpr_index = (1:K) + (K*D + choose(K,2))
bs_tpr_index = (1:K) + (K*(D+1) + choose(K,2))
bs_fpr_index = (1:K) + (K*(D+2) + choose(K,2))

log_JointDist = function(THETA){
      Prob_L.unique = XbetaTologPrL(K, X.unique, THETA[Beta_index], THETA[theta2_index], 
            LUmat, MuMat, theta1.init=NULL, normalize = TRUE)
      
      ProbMat.M_L = exp(log_ProbMat_MSSBSgivenL(K, Lmat.withZero, MSS, MBS, 
            THETA[ss_tpr_index], THETA[bs_tpr_index], THETA[bs_fpr_index]))
      ProbMat.L = Prob_L.unique[X.index,]

      Prob_M = as.vector(apply(ProbMat.M_L * ProbMat.L, 1, sum))

      log_Prob_case = sum(log(Prob_M), na.rm=TRUE)
      log_Prob = log_Prob_case + log_Prob_ctrl(K, N.ctrl, MBS.ctrl, THETA[bs_fpr_index]) + 
      Prior.BetaTheta2(THETA[Beta_index], THETA[theta2_index], 
            varbeta=varbeta, vartheta=vartheta, mutheta=mutheta, pmix=pmix, logscale=TRUE) +
      sum(dbeta(x=THETA[-(1:(K*D + choose(K,2)))], shape1=c(aa,cc,ee), shape2=c(bb,dd,ff), log=TRUE))
      return(-log_Prob)      
}

log_JointDist(THETA)


M.result = optim(par=THETA, fn=log_JointDist, method="Nelder", control=list(trace=2, maxit=2000))
cbind(THETA, round(M.result$par,3))


tmp =DirectMAP(K=5, Smax=5, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)

