library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
library(nleqslv)
library(BB)
library(mgcv)

sourceCpp("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/AQE_Components.cpp")
#---------------------------------------------------------------------------
logit = function(x)
{
      log(x/(1-x))
}

inv.logit = function(x)
{
      1/(1+exp(-x))
}

deriv.logit = function(x)
{
      1/x + 1/(1-x)
}


########################################################################################
################################### Design Matrix Preparation ###################################
########################################################################################

AQE.DesignMatrix = function(K, Smax=3){
      index = cumsum(sapply(1:Smax, function(x) choose(K,x)))
      Lmat = matrix(NA, ncol=K, nrow=index[Smax])
      row1 = 1
      for (s in 1:Smax){
            Lmat[row1:index[s],] = t(combn(1:K,s, function(x){a=rep(0,K);a[x]=1;a}))
            row1 = index[s]+1
      }

      Umat = t(apply(Lmat, 1, function(l) apply(combn(l,2),2,prod)))

      MuMat = matrix(NA,nrow=K,ncol=index[Smax])
      for (i in 1:K)
      {
            MuMat[i,] = as.numeric(Lmat[,i]==1)
      }
      return(list(Lmat=Lmat, Umat=Umat, MuMat=MuMat, J1=index[Smax]))
}
# AQE.DesignMatrix(5,3)
# dim(AQE.DesignMatrix(30,4))


# Get then row index in Lmat given M.GS
getCellLabel = function(K, Lmatrix)
{
      lmat = rbind(rep(0,K), Lmatrix)
      label = apply(lmat,1,function(x) paste(as.character(x),collapse=""))
      label
}
# cellLabel = getCellLabel(3,dmat$Lmat)

getLindex = function(Li, label)
{
      mgs = paste(as.character(Li),collapse="")
      which(label==mgs)
}
# getLindex(c(0,1,0), cellLabel)


########################################################################################
################################### Reparameterization #################################
########################################################################################

XbetaTologPrL = function(K, X.unique, Beta, theta2, lumat, normalize = TRUE){
      X.index = attr(X.unique, "index")
      Beta = matrix(Beta, ncol=K)
      theta1.unique = X.unique%*%Beta

      if (normalize == TRUE) {
            cellprobs.unique = t(apply(theta1.unique, 1, function(x) {c(1, exp(lumat%*%c(x, theta2)))/(1+sum(exp(lumat%*%c(x, theta2))))} ) )
      } else if (normalize == FALSE) {
            cellprobs.unique = t(apply(theta1.unique, 1, function(x) {c(1, exp(lumat%*%c(x, theta2))) } ) )
      } else {
             cellprobs.unique = t(apply(theta1.unique, 1, function(x) {c(0, lumat%*%c(x, theta2)) } ) )
      }
      return(cellprobs.unique)
}
# XbetaTologPrL(K=3, X.unique, Beta, c(-1,-0.5, 0.8), lumat, normalize=FALSE)


########################################################################################
################################### Simulation Tools ###################################
########################################################################################

AQE.simulate = function(K=3, Smax=3, P=c(0.5,0.2), N=50, Beta=NULL, theta2=c(-1,-0.5, 0.8), 
                                    ss_tpr, bs_tpr, bs_fpr, seed = 123, tol=1e-8) 
{
      set.seed(seed)

      D = length(P) + 1
      dmat = AQE.DesignMatrix(K,Smax)
      lumat = cbind(dmat$Lmat, dmat$Umat)
      lmat.withZero = rbind(rep(0,K), dmat$Lmat)

      if (length(Beta)<1 & D == 3) {
            Beta = matrix(NA, nrow=D, ncol=K)
            mu0 = c(0.3,0.45,0.2, 0.7, 0.2)
            Beta[1,] = logit(mu0[1:K])
            Beta[2,] = rbinom(K,1,0.5)-0.5
            Beta[3,] = rbinom(K,1,0.4)*0.5
      }

      X = matrix(NA, nrow=N, ncol=D)
      X[,1] = 1
      for (p in 2:D) {
            X[,p] = rbinom(N, 1, P[p-1])
      }
      X.unique = uniquecombs(X)
      X.index = attr(X.unique, "index")

      theta1.unique = X.unique%*%Beta
      cellprobs.unique = t(apply(theta1.unique, 1, function(x) {c(1, exp(lumat%*%c(x, theta2)))/(1+sum(exp(lumat%*%c(x, theta2))))} ) )

      
      Mu.unique = t(apply(cellprobs.unique, 1, function(x) {dmat$MuMat%*%x[-1]} ) ) 
      
      L = t(sapply(X.index, function(x) {t(rmultinom(1,1,cellprobs.unique[x,]))%*%lmat.withZero } ) )

      LtoM = function(L,TPR,FPR){
            k = ncol(L)
            notL = 1-L
            #registerDoMC(detectCores()) 
            M= foreach(i=1:nrow(L), .combine=rbind) %dopar% {
                  rbinom(k,1,L[i,]*TPR+notL[i,]*FPR)
            }
            M
      }

      # ss_tpr = c(0.05, 0.12, 0.08); bs_tpr = c(0.8,0.6, 0.7); bs_fpr = c(0.5, 0.55, 0.4)
      MSS.case = LtoM(L, ss_tpr, 0)
      MBS.case = LtoM(L, bs_tpr, bs_fpr)
      MBS.ctrl = t(matrix(rbinom(N*K, 1, bs_fpr), nrow=K))

      return(list(X=X, Beta=Beta, L=L, MSS.case=MSS.case, MBS.case=MBS.case, MBS.ctrl=MBS.ctrl, Mu.unique = Mu.unique, cellprobs=cellprobs.unique))
}
# dmat = AQE.DesignMatrix(3,3); lumat = cbind(dmat$Lmat, dmat$Umat); cellLabel = getCellLabel(3,dmat$Lmat)
# ss_tpr = c(0.05, 0.12, 0.08); bs_tpr = c(0.8,0.6, 0.7); bs_fpr = c(0.5, 0.55, 0.4)
# Data = AQE.simulate(N=2000, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=12345)



########################################################################################
################## Hamiltonian MC within Gibbs with Data Augmentation ###################
########################################################################################
Prior.BetaTheta2 = function(Beta, Theta2, varbeta, vartheta, mutheta, pmix, logscale=TRUE){
      pBeta = dnorm(as.vector(Beta), 0, varbeta, logscale)
      pTheta.a = dnorm(Theta2, mutheta, vartheta, FALSE)
      pTheta.b = dnorm(Theta2, 0, vartheta, FALSE)
      if (logscale) {
            priorvalue = sum(pBeta) + sum(log(pmix*pTheta.a + (1-pmix)*pTheta.b))
      } else {
            priorvalue = prod(pBeta)*prod(pmix*pTheta.a + (1-pmix)*pTheta.b)
      }
      return(priorvalue)
}
# Prior.BetaTheta2(Beta, Theta2, varbeta=20, vartheta=2, mutheta=-1, pmix=0.7, logscale=TRUE)


next_PositiveRates = function(K, L, MSS, MBS.case, MBS.ctrl, a, b, c, d, e, f){
      A = B = C = D = E = F = rep(NA, K)
      Ncase = nrow(L)
      Nctrl = nrow(MBS.ctrl)

      for (k in 1:K){
            sumLSS = crossprod(L[,k], MSS[,k])
            sumL = sum(L[,k])
            A[k] = sumLSS + a[k]
            B[k] = sumL - sumLSS + b[k]

            sumLBS = crossprod(L[,k], MBS.case[,k])
            C[k] = sumLBS + c[k]
            D[k] = sumL - sumLBS + d[k]

            sumBS = sum(MBS.case[,k])
            sumBSct = sum(MBS.ctrl[,k])
            E[k] = sumBS - sumLBS + e[k] + sumBSct
            F[k] = Ncase - sumBS - sumL + sumLBS + Nctrl - sumBSct 
      }
      ss_tpr.new = rbeta(K, A, B)
      bs_tpr.new = rbeta(K, C, D)
      bs_fpr.new = rbeta(K, E, F)
      return(rbind(ss_tpr.new, bs_tpr.new, bs_fpr.new))
}
# next_PositiveRates(K, Data$L, Data$MSS.case, Data$MBS.case, Data$MBS.ctrl, rep(1,3), rep(1,3),rep(1,3),rep(1,3),rep(1,3),rep(1,3))


## HMC step for beta & theta2
next_Lindex = function(K, Lmat.withZero, MSS, MBS.case, ss_tpr, bs_tpr, bs_fpr, X.index, LU_Theta.unique){
      log_ProbMat.a = log_ProbMat_MSSBSgivenL(K, Lmat.withZero, MSS, MBS.case, ss_tpr, bs_tpr, bs_fpr)
      log_ProbMat.b = LU_Theta.unique[X.index,]
      ProbMat = exp(log_ProbMat.a + log_ProbMat.b)
      indexSet = 1:nrow(Lmat.withZero)
      Lindex = apply(ProbMat,1, function(x) indexSet%*%rmultinom(1,1,x) )
      return(Lindex)
}
# log_ProbMat.a = log_ProbMat_MSSBSgivenL(K, Lmat.withZero, Data$MSS.case, Data$MBS.case, ss_tpr, bs_tpr, bs_fpr)
# X.index=attr(uniquecombs(Data$X), "index")
# LU_Theta.unique = XbetaTologPrL(3, X.unique, Data$Beta, theta2, lumat, dmat$MuMat, theta1.init=NULL, normalize = "linear")
# Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
# next_Lindex(K, Lmat.withZero, Data$MSS.case, Data$MBS.case, ss_tpr, bs_tpr, bs_fpr, X.index, LU_Theta.unique)


Hamiltonian_U = function(q, K, D, L.index, LUmat, X.unique, hyperPars){
      Beta = matrix(q[1:(K*D)], nrow=D, ncol=K)
      Theta2 = q[-(1:(K*D))]
      X.index = attr(X.unique, "index")

      U1.mat = XbetaTologPrL(K, X.unique, Beta, Theta2, LUmat, normalize = "linear")
      NormalizingConst = as.vector(apply(U1.mat, 1, function(x) {sum(exp(x)) } ) )

      data.index = cbind(L.index, X.index)
      U1 = sum(apply(data.index, 1, function(x) {U1.mat[x[2], x[1]] - log(NormalizingConst[x[2]]) })) 
      U2 = Prior.BetaTheta2(Beta, Theta2, hyperPars$varbeta, hyperPars$vartheta, hyperPars$mutheta, hyperPars$pmix, TRUE)
      U = -U1-U2
      return(U)
}
# Lindex = apply(Data$L,1,function(x) {getLindex(x, cellLabel)} )
# X.unique = uniquecombs(Data$X)
# theta2=c(-1,-0.5, 0.8)
# q = c(as.vector(Data$Beta), theta2)
# U = function(q) { Hamiltonian_U(q, K=3, D=3, Lindex, lumat, X.unique, hyperPars=data.frame(varbeta=20, vartheta=10, mutheta=-0.1, pmix=1)) }
# attr( numericDeriv(quote(U(q)), c("q")) , "gradient" ) -> g
# optim(par=q, fn=U, method="BFGS", control=list(trace=2))


HMC = function (U, epsilon, Nstep, current_q) {
      q = current_q
      p = rnorm(length(q),0,1) # independent standard normal variates 
      current_p = p
      # Make a half step for momentum at the beginning 
      grad_Uq = attr( numericDeriv(quote(U(q)), c("q")) , "gradient" )
      p = p - epsilon * grad_Uq/2
      # Alternate full steps for position and momentum
      for (i in 1:Nstep){
            # Make a full step for the position
            q = q + epsilon * p
            # Make a full step for the momentum, except at end of trajectory
            grad_Uq = tryCatch(attr( numericDeriv(quote(U(q)), c("q")) , "gradient" ) , error = function(c) stop("error finding gradient of U."))
            if (i!=Nstep) p = p - epsilon * grad_Uq
      }
      # Make a half step for momentum at the end.
      grad_Uq = tryCatch(attr( numericDeriv(quote(U(q)), c("q")) , "gradient" ) , error = function(c) stop("error finding gradient of U."))
      p = p - epsilon * grad_Uq/2
      # Negate momentum at end of trajectory to make the proposal symmetric
      p = -p
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_U = U(current_q) 
      current_K = sum(current_p^2) / 2 
      proposed_U = U(q)
      proposed_K = sum(p^2) / 2
      # Accept or reject the state at end of trajectory, returning either 
      # the position at the end of the trajectory or the initial position
      Hraio = exp(current_U-proposed_U+current_K-proposed_K)
      #print(paste0("Hratio: ", Hraio)); print(paste0(current_q, " vs. ", round(q,3))) ; 
      print(paste0(round(current_U,4), " vs. ", round(proposed_U,4)))
      if (runif(1) < Hraio){
            return (c(1,q)) # accept 
      } else {
            return (c(0,current_q)) # reject 
      }
}
# HMC(U, 0.05, 5, q)

# dmat = AQE.DesignMatrix(3,3); lumat = cbind(dmat$Lmat, dmat$Umat); cellLabel = getCellLabel(3,dmat$Lmat)
# ss_tpr = c(0.05, 0.12, 0.08); bs_tpr = c(0.8,0.6, 0.7); bs_fpr = c(0.5, 0.55, 0.4)
# Data = AQE.simulate(N=100, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=12345)
# Lindex = apply(Data$L,1,function(x) {getLindex(x, cellLabel)} )
# X.unique = uniquecombs(Data$X)
# theta2=c(-1,-0.5, 0.8)
# q = c(as.vector(Data$Beta), theta2)
# U = function(q) { Hamiltonian_U(q, K=3, D=3, Lindex, lumat, dmat$MuMat, X.unique, hyperPars=data.frame(varbeta=20, vartheta=10, mutheta=-0.1, pmix=0.7)) }
# attr( numericDeriv(quote(U(q)), c("q")) , "gradient" ) -> g
# optim(par=q, fn=U, method="BFGS")



# -------------------------------------------------------------------------------------------------------------------------- #

GetPosterior = function(K, Smax, HasGS=FALSE, 
                        DataList, # MSS, MBS, MBS.ctrl, N.ctrl, 
                        ParInit, #  beta, theta2, ss_tpr, bs_tpr, bs_fpr, 
                        HyperPar, # BetaPar, ThetaPar, a to f
                        Control ) # niter, nburn, epsilon, Nstep)
{
      # data & parameter input
      dmat = AQE.DesignMatrix(K,Smax)
      LUmat = cbind(dmat$Lmat, dmat$Umat)
      MuMat = dmat$MuMat
      J1 = dmat$J1
      Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
      cellLabel = apply(Lmat.withZero,1,function(x) paste(as.character(x),collapse=""))

      if (HasGS) {
            I_GS = DataList$I_GS
            MGS = DataList$MGS           
      }
      MSS = DataList$MSS.case
      MBS = DataList$MBS.case
      MBS.ctrl = DataList$MBS.ctrl
      X = DataList$X
      X.unique = uniquecombs(X)
      X.index = attr(X.unique, "index")

      D = ncol(X)
      N.case = nrow(MBS)
      N.ctrl = nrow(MBS.ctrl)

      iter = Control$iter
      burnin = Control$burnin
      epsilon = Control$epsilon
      Nstep = Control$Nstep

      aa = HyperPar$a; bb = HyperPar$b; cc = HyperPar$c;
      dd = HyperPar$d; ee = HyperPar$e; ff = HyperPar$f;
      BetaThetaPar = data.frame(varbeta=HyperPar$varbeta, vartheta=HyperPar$vartheta, 
                              mutheta=HyperPar$mutheta, pmix=HyperPar$pmix)

      # Initialize
      beta.posterior = matrix(NA, nrow=iter+burnin, ncol=K*D)
      theta2.posterior = matrix(NA, nrow=iter+burnin, ncol=choose(K,2))
      ss_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_fpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      Lindex.chain = matrix(NA, nrow=iter+burnin, ncol=N.case)

      beta.posterior[1,] = ParInit$beta
      theta2.posterior[1,] = ParInit$theta2
      ss_tpr.posterior[1,] = ParInit$ss_tpr
      bs_tpr.posterior[1,] = ParInit$bs_tpr
      bs_fpr.posterior[1,] = ParInit$bs_fpr


      LU_Theta.unique = XbetaTologPrL(K, X.unique, beta.posterior[1, ], theta2.posterior[1, ], LUmat, normalize = "linear")
      Lindex.chain[1, ] = next_Lindex(K, Lmat.withZero, MSS, MBS, ss_tpr.posterior[1,], 
                                    bs_tpr.posterior[1,], bs_fpr.posterior[1,], X.index, LU_Theta.unique)

      accept_track = numeric(iter+burnin)
      accept_track[1] = 0
      #loglik_track = numeric(iter+burnin)
      
      message("Start Hamiltonian MC within Gibbs sampling...")

      for (i in 2:(iter+burnin))
      {
            # TPRs, FPRs
            parm.a = next_PositiveRates(K, Lmat.withZero[Lindex.chain[i-1,],], 
                                          MSS, MBS, MBS.ctrl, aa, bb, cc, dd, ee, ff)
            ss_tpr.posterior[i,] = parm.a[1,]
            bs_tpr.posterior[i,] = parm.a[2,]
            bs_fpr.posterior[i,] = parm.a[3,]

            # (beta, theta2)
            current_q = c(beta.posterior[i-1, ], theta2.posterior[i-1, ])
            U = function(q) { Hamiltonian_U(q, K, D, Lindex.chain[i-1,], LUmat, X.unique, hyperPars=BetaThetaPar) }
            parm.b = tryCatch( HMC(U, epsilon, Nstep, current_q), error = function(c) "Incompatible Hamiltonian Proposal. Reject." )
            if (is.character(parm.b)){
                  print(paste0(i, ": ", parm.b))
                  accept_track[i] = 0
                  beta.posterior[i, ] = beta.posterior[i-1,]
                  theta2.posterior[i, ] = theta2.posterior[i-1,]
            } else {
                  accept_track[i] = parm.b[1]
                  beta.posterior[i, ] = parm.b[2:(K*D+1)]
                  theta2.posterior[i, ] = parm.b[-(1:(K*D+1))]
            }
            #print(c(i, accept_track[i]))

            # L
            if (accept_track[i]>0) {
                  LU_Theta.unique = XbetaTologPrL(K, X.unique, beta.posterior[i, ], theta2.posterior[i, ], LUmat, normalize = "linear")
            }
            Lindex.chain[i, ] = next_Lindex(K, Lmat.withZero, MSS, MBS, ss_tpr.posterior[i,], 
                                    bs_tpr.posterior[i,], bs_fpr.posterior[i,], X.index, LU_Theta.unique)

            # messages
            if (i%%25 == 0) 
            {
                  print(round(c(ss_tpr.posterior[i,], bs_tpr.posterior[i,], bs_fpr.posterior[i,]),3))
                  print(round(c(beta.posterior[i,]),3))
                  print("Trailing 200 sample mean:")
                  print(round(apply(beta.posterior[max(2, i-200):i,],2,mean),3))
            }
            if (i%%10 == 0) 
            {
                  print("Average accept rate:")
                  print(round(mean(accept_track[2:i]),3) )
                  #print(apply(accept_track[1:i,],2,mean))
                  print("Trailing 20 accept rate:")
                  print(round(mean(accept_track[max(2, i-19):i]),3) )
                  #print(apply(accept_track[max(1, i-19):i,],2,mean))
            }
      }
      result = list(beta.posterior = beta.posterior[-(1:burnin),], theta2.posterior = theta2.posterior[-(1:burnin),],
            ss_tpr.posterior = ss_tpr.posterior[-(1:burnin),], bs_tpr.posterior = bs_tpr.posterior[-(1:burnin),], 
            bs_fpr.posterior = bs_fpr.posterior[-(1:burnin),], accept_track = accept_track[-(1:burnin)], 
            Lindex = Lindex.chain[-(1:burnin),] )
      return(result)
}

# Data = AQE.simulate(N=500, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr, seed=12345)
# PAR0 = list(beta=as.vector(Data$Beta), theta2 = theta2, ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr)
# HYPAR = list(a=rep(1,K), b=rep(10,K), c=rep(2,K), d=rep(1,K), e=rep(1,K), f=rep(1,K), 
#             varbeta=20, vartheta=15, mutheta=-0.1, pmix=1)
# ctrlpar = data.frame(burnin=1, iter=1000, epsilon=0.03, Nstep=5)
# tmp = GetPosterior(K=3, Smax=3, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)

# round(as.vector(Data$Beta),3)
# round(apply(tmp$beta.posterior[-(1:500),],2,mean),3)
# theta2
# round(apply(tmp$theta2.posterior[-(1:500),],2,mean),3)
# round(apply(tmp$ss_tpr.posterior[-(1:500),],2,mean),3)
# bs_tpr
# round(apply(tmp$bs_tpr.posterior[-(1:500),],2,mean),3)

# table(tmp$Lindex[-(1:5),180]); getLindex(Data$L[180,], cellLabel)


########################################################################################
#################################### EM solution  #####################################
########################################################################################

# XbetaTologPrL(K, X.unique, beta.posterior[1, ], theta2.posterior[1, ], LUmat, MuMat, theta1.init=NULL, normalize = "linear")

GetEMsolution = function(K, Smax, HasGS=FALSE, 
                        DataList, # MSS, MBS, MBS.ctrl, N.ctrl, 
                        ParInit, #  beta, theta2, ss_tpr, bs_tpr, bs_fpr, 
                        HyperPar, # BetaPar, ThetaPar, a to f
                        Control ) # maxiter, tol
{
      # data & parameter input
      dmat = AQE.DesignMatrix(K,Smax)
      LUmat = cbind(dmat$Lmat, dmat$Umat)
      MuMat = dmat$MuMat
      J1 = dmat$J1
      Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
      cellLabel = apply(Lmat.withZero,1,function(x) paste(as.character(x),collapse=""))

      if (HasGS) {
            I_GS = DataList$I_GS
            MGS = DataList$MGS           
      }
      MSS = DataList$MSS.case
      MBS = DataList$MBS.case
      MBS.ctrl = DataList$MBS.ctrl
      X = DataList$X
      X.unique = uniquecombs(X)
      X.index = attr(X.unique, "index")

      D = ncol(X)
      N.case = nrow(MBS)
      N.ctrl = nrow(MBS.ctrl)

      maxiter = Control$maxiter
      tol = Control$tol

      aa = HyperPar$a; bb = HyperPar$b; cc = HyperPar$c;
      dd = HyperPar$d; ee= HyperPar$e; ff = HyperPar$f;
      varbeta=HyperPar$varbeta; vartheta=HyperPar$vartheta 
      mutheta=HyperPar$mutheta; pmix=HyperPar$pmix

      # Initialize
      THETA.track = matrix(NA, nrow=maxiter, ncol=K*D + choose(K,2) + 3*K)
      THETA.track[1,] = c(ParInit$beta, ParInit$theta2, ParInit$ss_tpr, ParInit$bs_tpr, ParInit$bs_fpr)
      Beta_index = 1:(K*D)
      theta2_index = (1:choose(K,2)) + (K*D)
      ss_tpr_index = (1:K) + (K*D + choose(K,2))
      bs_tpr_index = (1:K) + (K*(D+1) + choose(K,2))
      bs_fpr_index = (1:K) + (K*(D+2) + choose(K,2))

      GetQfunc = function(THETA, THETA_old){
            LU_Theta.unique_old = XbetaTologPrL(K, X.unique, THETA_old[Beta_index], THETA_old[theta2_index], 
                                    LUmat, normalize = "linear")
            log_ProbMat.a = log_ProbMat_MSSBSgivenL(K, Lmat.withZero, MSS, MBS, 
                              THETA_old[ss_tpr_index], THETA_old[bs_tpr_index], THETA_old[bs_fpr_index])
            log_ProbMat.b = LU_Theta.unique_old[X.index,]
            ProbMat = exp(log_ProbMat.a + log_ProbMat.b)
            NormalizingConst_old = as.vector(apply(ProbMat, 1, sum))
            ProbMat_old = ProbMat/NormalizingConst_old

            LU_Theta.unique = XbetaTologPrL(K, X.unique, THETA[Beta_index], THETA[theta2_index], 
                                    LUmat, normalize = "linear")
            NormalizingConst = as.vector(apply(exp(LU_Theta.unique), 1, sum))
            log_ProbMat.a = log_ProbMat_MSSBSgivenL(K, Lmat.withZero, MSS, MBS, 
                              THETA[ss_tpr_index], THETA[bs_tpr_index], THETA[bs_fpr_index])
            log_ProbMat.b = LU_Theta.unique[X.index,] - log(NormalizingConst[X.index])
            log_ProbMat = log_ProbMat.a + log_ProbMat.b

            q = sum(log_ProbMat*ProbMat_old, na.rm=TRUE)
            Qfunc = q + log_Prob_ctrl(K, N.ctrl, MBS.ctrl, THETA[bs_fpr_index]) + 
                  Prior.BetaTheta2(THETA[Beta_index], THETA[theta2_index], 
                  varbeta=varbeta, vartheta=vartheta, mutheta=mutheta, pmix=pmix, logscale=TRUE) +
                  sum(dbeta(x=THETA[-(1:(K*D + choose(K,2)))], shape1=c(aa,cc,ee), shape2=c(bb,dd,ff), log=TRUE))
            return(Qfunc)      
      }

      Qvalue_track = rep(NA, maxiter)
      Qvalue_track[1] = Inf


      log_JointDist = function(THETA){
            Prob_L.unique = XbetaTologPrL(K, X.unique, THETA[Beta_index], THETA[theta2_index], 
                  LUmat, normalize = TRUE)
      
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

      logProb_track = rep(NA, maxiter)
      logProb_track[1] = log_JointDist(THETA.track[1,])

      for (i in 2:maxiter){
            # E-step
            current_Q = function(x) {
                  -GetQfunc(THETA=x, THETA_old=THETA.track[i-1,])
            }
            # M-step
            print(paste0(i, ": Optimizing Q function..."))
            M.result = optim(par=THETA.track[i-1,], fn=current_Q, method="BFGS", control=list(trace=1))
            THETA.track[i,] = M.result$par
            Qvalue_track[i] = M.result$value
            logProb_track[i] = log_JointDist(THETA.track[i,])

            print(THETA.track[i, 1:(K*D)]); print(logProb_track[i])

            if (abs(logProb_track[i]-logProb_track[i-1]) < tol){
                  message("EM algorithm converged.")
                  result = list(THETA.final = THETA.track[i,], THETA.track = THETA.track[1:i,], 
                        Qvalue.track=Qvalue_track[1:i], logProb.track=logProb_track[1:i])
                  return(result)
            }
      }
      message("Maximum iteration is reached.")
      result = list(THETA.final = THETA.track[i,], THETA.track = THETA.track[1:i,], Qvalue.track=Qvalue_track[1:i],
            logProb.track=logProb_track[1:i])
      return(result)
}

# ctrlpar = data.frame(maxiter=100, tol=1e-5)
# tmp2 = GetEMsolution(K=3, Smax=3, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)
log_JointDist = function(THETA){
      Prob_L.unique = XbetaTologPrL(K, X.unique, THETA[Beta_index], THETA[theta2_index], 
            LUmat, normalize = TRUE)
      
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



########################################################################################
################################## Newton solution  #####################################
########################################################################################

DirectMAP = function(K, Smax, HasGS=FALSE, 
                        DataList, # MSS, MBS, MBS.ctrl, N.ctrl, 
                        ParInit, #  beta, theta2, ss_tpr, bs_tpr, bs_fpr, 
                        HyperPar, # BetaPar, ThetaPar, a to f
                        Control ) # maxiter, tol
{
      # data & parameter input
      dmat = AQE.DesignMatrix(K,Smax)
      LUmat = cbind(dmat$Lmat, dmat$Umat)
      MuMat = dmat$MuMat
      J1 = dmat$J1
      Lmat.withZero = rbind(rep(0,K), dmat$Lmat)
      cellLabel = apply(Lmat.withZero,1,function(x) paste(as.character(x),collapse=""))

      if (HasGS) {
            I_GS = DataList$I_GS
            MGS = DataList$MGS           
      }
      MSS = DataList$MSS.case
      MBS = DataList$MBS.case
      MBS.ctrl = DataList$MBS.ctrl
      X = DataList$X
      X.unique = uniquecombs(X)
      X.index = attr(X.unique, "index")

      D = ncol(X)
      N.case = nrow(MBS)
      N.ctrl = nrow(MBS.ctrl)

      maxiter = Control$maxiter
      tol = Control$tol

      aa = HyperPar$a; bb = HyperPar$b; cc = HyperPar$c;
      dd = HyperPar$d; ee= HyperPar$e; ff = HyperPar$f;
      varbeta=HyperPar$varbeta; vartheta=HyperPar$vartheta 
      mutheta=HyperPar$mutheta; pmix=HyperPar$pmix

      # Initialize
      THETA.track = matrix(NA, nrow=maxiter, ncol=K*D + choose(K,2) + 3*K)
      THETA.track[1,] = c(ParInit$beta, ParInit$theta2, ParInit$ss_tpr, ParInit$bs_tpr, ParInit$bs_fpr)
      Beta_index = 1:(K*D)
      theta2_index = (1:choose(K,2)) + (K*D)
      ss_tpr_index = (1:K) + (K*D + choose(K,2))
      bs_tpr_index = (1:K) + (K*(D+1) + choose(K,2))
      bs_fpr_index = (1:K) + (K*(D+2) + choose(K,2))

      log_JointDist = function(THETA){
            Prob_L.unique = XbetaTologPrL(K, X.unique, THETA[Beta_index], THETA[theta2_index], 
                                    LUmat, normalize = TRUE)
            
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

      M.result = optim(par=THETA.track[1,], fn=log_JointDist, method="BFGS", control=list(trace=1))
      return(M.result)
}

# tmp3 =DirectMAP(K=3, Smax=3, DataList=Data, ParInit=PAR0, HyperPar=HYPAR, Control=ctrlpar)

