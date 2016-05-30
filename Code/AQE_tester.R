library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
library(nleqslv)
library(BB)
library(mgcv)

sourceCpp("/Users/dengdetian0603/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/AQE_Components.cpp")
source("/Users/dengdetian0603/Documents/JHSPH/Research/S.Zeger/LogLinearModel/Code/AQE_models_marginal.R")
#register()

#---------------------------------------------------------------------------
K=5; Smax=4
dmat = AQE.DesignMatrix(K,Smax)
cellLabel = getCellLabel(K,dmat$Lmat)

ss_tpr = c(0.05, 0.12, 0.08, 0.15, 0.1); bs_tpr = c(0.8,0.6,0.7,0.7,0.5); bs_fpr = c(0.5, 0.55, 0.4, 0.35, 0.45)
#theta2 = rep(0,10)
theta2=c(-0.1,-0.3, -0.4, -0.8, -0.45, -0.3, -0.4, -0.7, 0.4, -0.9)

Data = AQE.simulate(K=5, Smax, N=5000, theta2=theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=rnorm(1,1000,100))

PAR0 = list(beta=as.vector(Data$Beta), theta2 = theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr)
HYPAR = list(a=rep(1,K), b=rep(10,K), c=rep(2,K), d=rep(1,K), e=rep(1,K), f=rep(1,K), 
             varbeta=20, vartheta=15, mutheta=-0.1, pmix=0.7)
# ----------------------------------------------------------------------------

LUmat = cbind(dmat$Lmat, dmat$Umat)
MuMat = dmat$MuMat
Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
#cellLabel = apply(Lmat.withZero,1,function(x) paste(as.character(x),collapse=""))

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

# ----------------------------- debugging Cpp functions ------------------------------------------
# with Cpp functions, X.index must starts from 0
W = EM_GetWeights(K, Lmat.withZero, MSS, MBS, ss_tpr, bs_tpr, bs_fpr, X.index-1, X.unique, LUmat, 
              as.vector(Data$Beta), theta2)

rates = EM_UpdateRates(K, nrow(LUmat), N.case, N.ctrl, MSS, MBS, MBS.ctrl, W, X.index-1, Lmat.withZero,
               aa, bb, cc, dd, ee, ff)

new_par = EM_UpdateBetaTheta2(W, X.index-1, X.unique, LUmat, K, nrow(LUmat), D,
                    varbeta, vartheta, mutheta, c(as.vector(Data$Beta), theta2))

stepaway = rnorm(25, 0, 0.05)
old_par = c(as.vector(Data$Beta), theta2) +  stepaway
qfunc(W, X.index-1, X.unique, LUmat, K, nrow(LUmat), D, varbeta, vartheta, mutheta, old_par)

tmp = XbetaToProbL_unique(K, nrow(LUmat), X.unique, old_par[1:(D*K)], old_par[-(1:(D*K))], LUmat, 3)


# ----------------------------- classic EM algorithm -------------------------------------------
# stepaway = rnorm(25, 0, 0.5)
# new_par = c(as.vector(Data$Beta), theta2) +  stepaway
# new_rates = rates

# new_par = rnorm(25, -1, 0.5)

# for (i in 1:250){
#       old_beta = new_par[1:(K*D)]
#       old_theta2 = new_par[-(1:(K*D))]
#       old_par = new_par
#       print(paste("iter:", i))
#       print("E-step")
#       W = EM_GetWeights(K, Lmat.withZero, MSS, MBS, new_rates[1,], new_rates[2,], new_rates[3,], 
#                   X.index-1, X.unique, LUmat, old_beta, old_theta2)
#       new_rates = EM_UpdateRates(K, nrow(LUmat), N.case, N.ctrl, MSS, MBS, MBS.ctrl, W, X.index-1, Lmat.withZero,
#                aa, bb, cc, dd, ee, ff)
#       print("M-step")
#       new_par = EM_UpdateBetaTheta2(W, X.index-1, X.unique, LUmat, K, nrow(LUmat), D,
#                     varbeta, vartheta, mutheta, old_par)
#       print( sqrt(sum((old_par-new_par)^2)) )
#       print( sqrt(sum((c(as.vector(Data$Beta), theta2)-new_par)^2)) )
# }

# ----------------------------- accelerated EM algorithms ---------------------------------------

library(turboEM)
registerDoMC(detectCores())
FixMap = function(x){
      EM_update(x, K, D, Lmat.withZero, MSS, MBS, MBS.ctrl, X.index, X.unique, 
                     LUmat, N.case, N.ctrl, aa, bb, cc, dd, ee, ff, varbeta, vartheta, mutheta)
}

par0 = rnorm(25, -1, 0.5)
tmp2 = turboem(par = c(ss_tpr, bs_tpr, bs_fpr, par0), 
      fixptfn = FixMap, objfn = ObjFunc, method = "squarem", parallel = TRUE) 
      #boundary, pconstr = NULL, project = NULL, parallel = FALSE, ...,
      #control.method = replicate(length(method),list()), control.run = list())
(tmp$pars - c(ss_tpr, bs_tpr, bs_fpr,  as.vector(Data$Beta), theta2))
cbind(par0[1:15], tmp2$pars[16:30], as.vector(Data$Beta))
cbind(tmp$pars[1:15], c(ss_tpr, bs_tpr, bs_fpr))

ObjFunc = function(x){
      -Evaluate_log_distn(x,X.index-1, X.unique, LUmat, K, nrow(LUmat), D, Lmat.withZero, 
                   MSS, MBS, MBS.ctrl, varbeta, vartheta, mutheta, 
                   aa, bb, cc, dd, ee, ff)
}
ObjFunc(tmp$par[1,])
ObjFunc(c(ss_tpr, bs_tpr, bs_fpr,  as.vector(Data$Beta), theta2))
