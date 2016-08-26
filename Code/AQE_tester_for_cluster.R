library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
library(nleqslv)
library(BB)
library(mgcv)
library(turboEM)


sourceCpp("/users/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_Components.cpp")
source("/users/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_models_marginal.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1){
      Smax = as.numeric(args[1])
      D = as.numeric(args[2])
      n_rep = as.numeric(args[3])
      exch = as.numeric(args[4])
} else {
      Smax = 3
      D = 1
      n_rep = 2
      exch = 1
}
K= 5 
n_each = 2000
maxCore = 20

print(paste0("Smax=",Smax," D=",D," exch=",exch," maxCore=",maxCore))
registerDoMC(min(detectCores(), maxCore))

#---------------------------------------------------------------------------
dmat = AQE.DesignMatrix(K,Smax)
cellLabel = getCellLabel(K,dmat$Lmat)

LUmat = cbind(dmat$Lmat, dmat$Umat)
MuMat = dmat$MuMat
Lmat.withZero = rbind(rep(0,K), dmat$Lmat)

ss_tpr = c(0.05, 0.12, 0.08, 0.15, 0.1); bs_tpr = c(0.8,0.6,0.7,0.7,0.5); bs_fpr = c(0.5, 0.55, 0.4, 0.35, 0.45)
if (exch>0.5) {
      theta2 = rep(-0.4, choose(K,2))
} else {
      theta2 = c(-0.1,-0.3, -0.4, -0.8, -0.45, -0.3, -0.4, -0.7, 0.4, -0.9)
}

theta2 = theta2[1:choose(K,2)]
ss_tpr = ss_tpr[1:K]
bs_tpr = bs_tpr[1:K]
bs_fpr = bs_fpr[1:K]


SEED = 1066.887
BETA0 = AQE.simulate(K=K, Smax=Smax, N=50, P=0.5, theta2=theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=SEED)$Beta
BETA0 = matrix( BETA0[1:D, 1:K], nrow=D )

print(BETA0)

PAR0 = list(beta=as.vector(BETA0), theta2 = theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr)
HYPAR = list(a=rep(1,K), b=rep(10,K), c=rep(2,K), d=rep(1,K), e=rep(1,K), f=rep(1,K), 
             varbeta=20, vartheta=5, mutheta=-0.3, pmix=0.7)
# ----------------------------------------------------------------------------

model_fits = foreach(k = 1:n_rep, .combine=rbind) %dopar% { 
      print(k)

      print("Simulating data...")
      Data = AQE.simulate(K=K, Smax=Smax, N=n_each, P=NULL, Beta=BETA0, theta2=theta2, 
                        ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=21205+7*k)

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

      # ----------------------------- accelerated EM algorithms ---------------------------------------

      FixMap = function(x){
            EM_update(x, K, D, Lmat.withZero, MSS, MBS, MBS.ctrl, X.index, X.unique, 
                           LUmat, N.case, N.ctrl, aa, bb, cc, dd, ee, ff, varbeta, vartheta, mutheta)
      }

      FixMap_exch = function(x){
            # EM_update_exch(x, K, D, Lmat.withZero, MSS, MBS, MBS.ctrl, X.index, X.unique, 
            #             LUmat, N.case, N.ctrl, aa, bb, cc, dd, ee, ff, varbeta, vartheta, mutheta)
            EM_update_exch_UpdateBetaWithBSonly(x, K, D, Lmat.withZero, MSS, MBS, MBS.ctrl, X.index, X.unique, 
                        LUmat, N.case, N.ctrl, aa, bb, cc, dd, ee, ff, varbeta, vartheta, mutheta)
      }

      ObjFunc = function(x){
            -Evaluate_log_distn(x,X.index-1, X.unique, LUmat, K, nrow(LUmat), D, Lmat.withZero, 
                         MSS, MBS, MBS.ctrl, varbeta, vartheta, mutheta, 
                         aa, bb, cc, dd, ee, ff)
      }

      ObjFunc_exch = function(x){
            -Evaluate_log_distn_exch(x,X.index-1, X.unique, LUmat, K, nrow(LUmat), D, Lmat.withZero, 
                   MSS, MBS, MBS.ctrl, varbeta, vartheta, mutheta, 
                   aa, bb, cc, dd, ee, ff)
      }

      ObjFunc_exch_BSonly = function(x){
            # TODO:
            -Evaluate_log_distn_exch(x,X.index-1, X.unique, LUmat, K, nrow(LUmat), D, Lmat.withZero, 
                   MSS, MBS, MBS.ctrl, varbeta, vartheta, mutheta, 
                   aa, bb, cc, dd, ee, ff)
      }

      par0 = c( as.vector(Data$Beta), theta2[1])
      par0 = par0 + rnorm(length(par0))
      
      print("EM optimizing...")
      print(paste0("Starting value: ", par0))
      time0 = proc.time()
      tmp = turboem(par = c(ss_tpr, bs_tpr, bs_fpr, par0), 
            fixptfn = FixMap_exch, #objfn = ObjFunc_exch,
            method = "squarem", 
            parallel = TRUE,
            control.run = list(
                  tol = 1e-6,
                  maxiter = 500,
                  trace=TRUE)) 
      print(proc.time()-time0)

      print(cbind(inv.logit(tmp$pars[(3*K+1):((3+D)*K)]),
                  inv.logit(as.vector(Data$Beta)),
                  inv.logit(tmp$pars[(3*K+1):((3+D)*K)]) - inv.logit(as.vector(Data$Beta))))
      print(cbind(tmp$pars[1:(3*K)],
                  c(ss_tpr, bs_tpr, bs_fpr),
                  tmp$pars[1:(3*K)]- c(ss_tpr, bs_tpr, bs_fpr)) )
      
      tmp$par
}

print(paste0("K5-Smax",Smax,"-D",D,"-exch",exch))

rates_fit = model_fits[,1:(3*K)]
beta_fit = model_fits[,(3*K+1):((3+D)*K)]
theta2_fit = model_fits[,-(1:((3+D)*K))]

beta_bias = inv.logit(matrix(apply(beta_fit, 2, mean), ncol=K)) - inv.logit(BETA0)
rates_bias = apply(rates_fit, 2, mean) - c(ss_tpr, bs_tpr, bs_fpr)
#theta2_bias = apply(theta2_fit, 2, mean) - PAR0$theta2

print("bias of beta:")
print(beta_bias)

save.image(file=paste0("K5-Smax",Smax,"-D",D,"-exch",exch,"-0808.Rdata"))
