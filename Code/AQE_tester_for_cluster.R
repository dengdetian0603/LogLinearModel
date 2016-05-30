library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
library(nleqslv)
library(BB)
library(mgcv)
library(turboEM)


sourceCpp("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_Components.cpp")
source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_models_marginal.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1){
      Smax = as.numeric(args[1])
      n_each = as.numeric(args[2])
      n_rep = as.numeric(args[3])
      maxCore = as.numeric(args[4])
} else {
      Smax = 4
      n_each = 1000
      n_rep = 50
      maxCore = 50
}

registerDoMC(min(detectCores(), maxCore))

#---------------------------------------------------------------------------
K=5; #Smax=4
dmat = AQE.DesignMatrix(K,Smax)
cellLabel = getCellLabel(K,dmat$Lmat)

LUmat = cbind(dmat$Lmat, dmat$Umat)
MuMat = dmat$MuMat
Lmat.withZero = rbind(rep(0,K), dmat$Lmat)

ss_tpr = c(0.05, 0.12, 0.08, 0.15, 0.1); bs_tpr = c(0.8,0.6,0.7,0.7,0.5); bs_fpr = c(0.5, 0.55, 0.4, 0.35, 0.45)
#theta2 = rep(0,10)
theta2=c(-0.1,-0.3, -0.4, -0.8, -0.45, -0.3, -0.4, -0.7, 0.4, -0.9)
SEED = 914.5116

BETA0 = AQE.simulate(K=5, Smax=Smax, N=5000, theta2=theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=SEED)$Beta
PAR0 = list(beta=as.vector(BETA0), theta2 = theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr)
HYPAR = list(a=rep(1,K), b=rep(10,K), c=rep(2,K), d=rep(1,K), e=rep(1,K), f=rep(1,K), 
             varbeta=20, vartheta=15, mutheta=-0.1, pmix=0.7)
# ----------------------------------------------------------------------------

model_fits = foreach(k = 1:n_rep, .combine=rbind) %dopar% { 
      print(k)

      print("Simulating data...")
      Data = AQE.simulate(K=5, Smax=Smax, N=n_each, Beta=BETA0, theta2=theta2, 
                        ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=21205+10*k)

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

      ObjFunc = function(x){
            -Evaluate_log_distn(x,X.index-1, X.unique, LUmat, K, nrow(LUmat), D, Lmat.withZero, 
                         MSS, MBS, MBS.ctrl, varbeta, vartheta, mutheta, 
                         aa, bb, cc, dd, ee, ff)
      }

      par0 = rnorm(25, -1, 0.5)
      
      print("EM optimizing...")
      time0 = proc.time()
      tmp = turboem(par = c(ss_tpr, bs_tpr, bs_fpr, par0), 
            fixptfn = FixMap, objfn = ObjFunc, method = "squarem", parallel = FALSE) 
            #boundary, pconstr = NULL, project = NULL, parallel = FALSE, ...,
            #control.method = replicate(length(method),list()), control.run = list())
      #print( tmp$pars - c(ss_tpr, bs_tpr, bs_fpr,  as.vector(Data$Beta), theta2) )
      #print( cbind(par0[1:15], tmp$pars[16:30], as.vector(Data$Beta)) )
      #print( cbind(tmp$pars[1:15], c(ss_tpr, bs_tpr, bs_fpr)) )
      print(proc.time()-time0)

      tmp$par
}

save(model_fits, BETA0, PAR0, HYPAR, n_each, Smax, file=paste("K5-Smax_",Smax,"-n",n_each,".Rdata", sep="" ) )
