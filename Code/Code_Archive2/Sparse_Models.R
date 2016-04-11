library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
# dependent on ToolBox_PERCH.R

GenMatrices.Sparse = function(K, Smax=3)
{
      L.all = as.matrix(AllComb(K))[-1,]
      S = apply(L.all,1,sum)
      index = order(S)
      L.all = L.all[index,]
      S = S[index]
      
      Sallowed = which(S<=Smax)
      L.all = L.all[Sallowed,]
      S = S[Sallowed]
      J1 = length(S)
      
      #-- expand all possible combinations
      #L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
      #      expandL.LL(L.all[i,])
      #}
      
      MuMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            MuMat[i,] = as.numeric(L.all[,i]==1)
      }
      
      PiMat = matrix(NA,nrow=Smax,ncol=J1)
      for (i in 1:Smax)
      {
            PiMat[i,] = as.numeric(S==i)
      }
      return(list(MuMat=MuMat, PiMat=PiMat, J1=J1, Lmatrix=L.all))
}

# Get then row index in Lmat given M.GS
getCellLabel = function(K, pars)
{
      J = pars$J1 + 1
      lmat = rbind(rep(0,K), pars$Lmatrix[,1:K])
      label = apply(lmat,1,function(x) paste(as.character(x),collapse=""))
      label
}
# cellLabel = getCellLabel(3,Lmat)

getGSindex = function(MGSi, label)
{
      mgs = paste(as.character(MGSi),collapse="")
      which(label==mgs)
}
# getGSindex(MGS=c(0,1,0), cellLabel)



########################### Likelihood components ####################################

### simulate phi given constraint defined by mu and pi
rPhiGivenPiMu.Sparse = function(K, Smax, mu, pis, burnin=30, n=10, method="mirror", pars)
{
      # TODO: number of iter and output should be adaptive according to K
      pis = pis[1:(Smax+1)]
      if (abs(sum(pis)-1)>1e-6)
      {
            print(sum(pis))
            message("pi vector is not Sparse.")
            return(NA)
      }
      Phis = tryCatch(xsample( E = rbind(pars$PiMat,pars$MuMat), F = c(pis[-1],mu)/pis[1], 
                      G = diag(rep(1,pars$J1)), H = rep(0,pars$J1), tol=1e-7, 
                      iter = n, burnin = burnin, type = method, test=FALSE), error = function(e) "Incompatible constraints.")
      if (is.character(Phis)) return(NA)
      else return( tryCatch(Phis$X, error=function(e){ print(e); print("pi:");print(pis);print("mu:"); print(mu); return(NA)}))
}
#rPhiGivenPiMu.Sparse(K=3,Smax=3,mu=TrueMu.exact, pis=TruePi.exact, n=50, method="mirror",pars=Lmat)

### Pr[L = l_j| mu, pi]
Prob.LgivenMuPi = function(K, Smax, mu, pis, n=500, sparse.par)
{
      phis = rPhiGivenPiMu.Sparse(K=K,Smax=Smax,mu=mu, pis=pis, burnin=30, n=n, method="mirror",pars=sparse.par)
      if (is.na(phis)) return(NA)
      else {
            cell.prob = colMeans(phis)*pis[1]
            return(c(pis[1], cell.prob))
      }
}
# Prob.LgivenMuPi(n=3000, K=3, Smax=3, mu=TrueMu.exact, pis=TruePi.exact, sparse.par=Lmat)->tmp


### Pr[MSS_i, MBS_i|l_j]
Prob.MSSBSigivenLj = function(K, MSSi, MBSi, Lj, ss_tpr, bs_tpr, bs_fpr, logscale=TRUE)
{
      logprobs = foreach(k = 1:K, .combine='+') %do% {
            #p1 = ((Lj[k]*ss_tpr[k]^Lj[k])^MSSi[k])*(1 - ss_tpr[k])^(Lj[k]*(1-MSSi[k]))
            #p2 = ((bs_tpr[k]^Lj[k] * bs_fpr[k]^(1-Lj[k]))^MBSi[k]) * ((1-bs_tpr[k])^Lj[k] * (1-bs_fpr[k])^(1-Lj[k]))^(1-MBSi[k])

            log.p1 = ifelse( Lj[k] + MSSi[k] < 0.5, 0, 
                         MSSi[k]*(log(Lj[k]) + Lj[k]*log(ss_tpr[k])) + Lj[k]*(1-MSSi[k])*log(1-ss_tpr[k]) )
            log.p2 = MBSi[k]*(Lj[k]*log(bs_tpr[k]) + (1-Lj[k])*log(bs_fpr[k])) + 
                        (1-MBSi[k])*(Lj[k]*log(1-bs_tpr[k]) + (1-Lj[k])*log(1-bs_fpr[k]))
            log.p1 + log.p2
      }
      if (logscale) return(logprobs)
      else return (exp(logprobs))
}
# t0=proc.time()
# Prob.MSSBSigivenLj(K=3, MSSi=M.SS[1,], MBSi=M.BS[1,], Lj=case.state[1,], ss_tpr=SS_TPR, bs_tpr=BS_TPR, bs_fpr=BS_FPR)
# log_Prob_MSSBSigivenLj(K=3, MSSi=M.SS[1,], MBSi=M.BS[1,], Lj=case.state[1,], ss_tpr=SS_TPR, bs_tpr=BS_TPR, bs_fpr=BS_FPR)
# proc.time()-t0



###################################### Likelihood #########################################

### Pr[MSS, MBS | mu, pi] Likelihood without GS measurements
Prob.SS_BS = function(K, MSS, MBS, I_GS, Pr_Lj, ss_tpr, bs_tpr, bs_fpr, sparse_pars, logscale=TRUE)
{
      J = sparse_pars$J1 + 1
      lmat = rbind(rep(0,K), sparse_pars$Lmatrix[,1:K])

      log.prob_i = foreach(i = which(I_GS<1), .combine='+') %dopar% {
            #TP.index = which(MSS[i,]>0)
            #nonZero.index = which(rowSums(lmat[,TP.index])==length(TP.index))

            #message(paste("+ i: ", i))
            #t0=proc.time()
            prob_ij = foreach(j = 1:J, .combine='+') %do% {
                  Pr_Lj[j] * exp(Prob.MSSBSigivenLj(K=K, MSSi=MSS[i,], MBSi=MBS[i,], Lj=lmat[j,], 
                                                ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr))
            }
            #message(paste("- i: ",i," time: ", (proc.time()-t0)[3]))
            log(prob_ij)
      }
      if (logscale) return(log.prob_i)
      else return(exp(log.prob_i))
}
# t0=proc.time()
# Prob.SS_BS(K=3, MSS=M.SS, MBS=M.BS, I_GS=I_GS, Pr_Lj=tmp, SS_TPR, BS_TPR, BS_FPR, sparse_pars=Lmat, logscale=TRUE)
# proc.time() - t0

# t0=proc.time()
# log_Prob_SS_BS(K=3, MSS=M.SS, MBS=M.BS, I_GS=I_GS, Pr_Lj=tmp, SS_TPR, BS_TPR, BS_FPR, Lall_matrix=rbind(rep(0,K), Lmat$Lmatrix[,1:K]))
# proc.time() - t0
# Prob.SS_BS(K=3, MSS=M.SS, MBS=M.BS, I_GS=I_GS, Pr_Lj=tmp, c(.99,.79,.90), c(.79,.89,.95), rep(0.01,3), sparse_pars=Lmat, logscale=TRUE)


### Pr[MGS, MSS, MBS | mu, pi] Likelihood with GS measurements
Prob.GS_SS_BS = function(K, MGS, MSS, MBS, I_GS, Pr_Lj, ss_tpr, bs_tpr, bs_fpr, logscale=TRUE, cell.label)
{
      log.prob_i = foreach(i = which(I_GS>0), .combine='+') %dopar% {
            log.prob_ik = foreach(k = 1:K, .combine='+') %do% {
                  p1 = ((MGS[i,k]*ss_tpr[k]^MGS[i,k])^MSS[i,k])*(1 - ss_tpr[k])^(MGS[i,k]*(1-MSS[i,k]))
                  p2 = ((bs_tpr[k]^MGS[i,k] * bs_fpr[k]^(1-MGS[i,k]))^MBS[i,k]) * ((1-bs_tpr[k])^MGS[i,k] * 
                        (1-bs_fpr[k])^(1-MGS[i,k]))^(1-MBS[i,k])
                  log(p1)+log(p2)
            }
            GSindx = match(paste(as.character(MGS[i,]),collapse=""), cell.label) ## this step might be slow for large J1
            p3 = Pr_Lj[GSindx]
            log.prob_ik+log(p3)
      }
      if (logscale) return(log.prob_i)
      else return(exp(log.prob_i))
}
# Prob.LgivenMuPi(n=6000, K=3, Smax=3, mu=TrueMu.exact, pis=TruePi.exact, sparse.par=Lmat)->tmp1
# Prob.LgivenMuPi(n=6000, K=3, Smax=3, mu=TrueMu.exact[c(1,3,2)], pis=TruePi.exact, sparse.par=Lmat)->tmp2
# Prob.GS_SS_BS(K=3, MGS=M.GS[downsample,], MSS=M.SS[downsample,], MBS=M.BS[downsample,], I_GS=I_GS[downsample], Pr_Lj=tmp1, SS_TPR, BS_TPR, BS_FPR, logscale=TRUE, cell.label)
# Prob.GS_SS_BS(K=3, MGS=M.GS, MSS=M.SS, MBS=M.BS, I_GS=I_GS, Pr_Lj=tmp, c(.8,.79,.80), c(.79,.69,.85), rep(0.3,3), logscale=TRUE)


### Likelihood for control data
Prob.ctrl = function(K, Nctrl, Mbs, bs_fpr, logscale=TRUE)
{
      log.prob = foreach(I = 1:(K*Nctrl), .combine='+') %dopar% {
            k = I%%K
            k = ifelse(k<1,K,k)
            i = (I-k)/K + 1
            log.prob_ik = Mbs[i,k]*log(bs_fpr[k]) + (1-Mbs[i,k])*log(1-bs_fpr[k])
            log.prob_ik
      }
      if (logscale) return(log.prob)
      else return(exp(log.prob))
}
# 
# t0=proc.time()
# Prob.ctrl(K=3, Nctrl=1000, Mbs=M.BS.ctrl, BS_FPR, TRUE)
# proc.time()-t0

# t0=proc.time()
# log_Prob_ctrl(K=3, Nctrl=1000, Mbs=M.BS.ctrl, BS_FPR)
# proc.time()-t0
# Prob.ctrl(K=3, Nctrl=1000, Mbs=M.BS.ctrl, c(0.1,0.1,0.1), TRUE)



###################################### priors #############################################
# partially informative prior on (mu, pi)
prior.MuPi.Sparse = function(mu, pis, Alpha, Smax, logscale=TRUE)
{
      K = length(mu)
      pis = pis[1:(Smax+1)]
      density.pi0toSmax = ddirichlet(x=pis, alpha=Alpha)
      q = 1 - pis[1]
      Q = crossprod(1:Smax, pis[-1])
      boundary1 = as.numeric(abs(sum(mu)-Q) < 1e-9)
      density.mu = ((1/q)^K)/dIrwinHall(x=Q, K=K, q=q)*prod(mu>0 & mu<q)*boundary1
      if (logscale)
      {
            return(log(density.mu) + log(density.pi0toSmax))
      }
      else
      {
            
            return(density.mu*density.pi0toSmax)
      }
}

# prior.MuPi.Sparse(mu=TrueMu.exact, pis=TruePi.exact, Alpha=c(40,10,20,10), Smax=3, logscale=TRUE)


# informative priors for TPRs and FPR
prior.TPR.FPR = function(ss_tpr, bs_tpr, bs_fpr, hyperPars, logscale=TRUE)
{# hyperPars needs to be a 3K x 2 matrix
      x = c(ss_tpr, bs_tpr, bs_fpr)
      alpha = hyperPars[,1]
      beta = hyperPars[,2]
      log.prior = sum(dbeta(x=x, alpha, beta, log=TRUE))
      if (logscale) return(log.prior)
      else return(exp(log.prior))
}

# hpar = matrix(c(rep(30,6),rep(2,3), c(1.6, 7.5, 3.5, 7.5, 3.5, 1.6, 198, 38, 18)), ncol=2)
# prior.TPR.FPR(ss_tpr=SS_TPR, bs_tpr=BS_TPR, bs_fpr=BS_FPR, hyperPars=hpar)
# prior.TPR.FPR(ss_tpr=c(0.7,0.8,0.8), bs_tpr=BS_TPR, bs_fpr=BS_FPR, hyperPars=hpar)


###################################### Full Joint density ##########################################
sourceCpp("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/ModelComponents.cpp")

Joint.Density.Cpp = function(K, Smax, HasGS=TRUE, MGS=NULL, I_GS=NULL, MSS, MBS, Mbs_ctrl, N_ctrl, 
                                    latent.cell.probs, mus, pis, ss_tpr, bs_tpr, bs_fpr, DirichPar, BetaPars, 
                                    logscale=TRUE, cell.labels=NULL, Lall_mat=NULL, verbatim=FALSE)
{

      lik.case.withGS = ifelse(HasGS, Prob.GS_SS_BS(K=K, MGS=MGS, MSS=MSS, MBS=MBS, I_GS=I_GS, 
                                    Pr_Lj=latent.cell.probs, ss_tpr, bs_tpr, bs_fpr, logscale=TRUE, cell.label=cell.labels),0)

      lik.case.noGS = log_Prob_SS_BS(K=K, MSS=MSS, MBS=MBS, I_GS=I_GS, 
                                    Pr_Lj=latent.cell.probs, ss_tpr, bs_tpr, bs_fpr, Lall_matrix=Lall_mat)

      lik.ctrl = log_Prob_ctrl(K=K, Nctrl=N_ctrl, Mbs=Mbs_ctrl, bs_fpr)
      
      prior.mupi = prior.MuPi.Sparse(mu=mus, pis=pis, Alpha=DirichPar, Smax=Smax, logscale=TRUE)
      prior.tfpr = prior.TPR.FPR(ss_tpr, bs_tpr, bs_fpr, hyperPars=BetaPars, logscale=TRUE)

      if (verbatim)
      {
      print(paste("GS Likelihood:",lik.case.withGS))
      print(paste("no GS likilihood: ", lik.case.noGS))
      print(paste("ctrl likilihood: ", lik.ctrl))
      print(paste("mu, pi priors: ", prior.mupi))
      print(paste("tfpr priors: ", prior.tfpr))
      }

      log.density = lik.case.withGS + lik.case.noGS + lik.ctrl + prior.mupi + prior.tfpr
      if (logscale) return(log.density)
      else return(exp(log.density))
}

############## Pre-generated parameters needed ###############
# Lmat = GenMatrices.Sparse(3,3)
# cellprob = Prob.LgivenMuPi(n=3000, K=3, Smax=3, mu=TrueMu.exact, pis=TruePi.exact, sparse.par=Lmat)
# hpar = matrix(c(rep(30,6),rep(2,3), c(1.6, 7.5, 3.5, 7.5, 3.5, 1.6, 198, 38, 18)), ncol=2)
# lall.mat = rbind(rep(0,K), Lmat$Lmatrix[,1:K])
# cell.label = apply(lall.mat,1,function(x) paste(as.character(x),collapse=""))
##############################################################
# Joint.Density.Cpp(K=3, Smax=3, MGS=M.GS[downsample,], I_GS=I_GS[downsample], MSS=M.SS[downsample,], MBS=M.BS[downsample,], Mbs_ctrl=M.BS.ctrl, N_ctrl=1000, 
#                    latent.cell.probs=cellprob, mus=TrueMu.exact, pis=TruePi.exact, ss_tpr=SS_TPR, bs_tpr=BS_TPR, bs_fpr=BS_FPR, 
#                    DirichPar=c(1,4,2,1), BetaPars=hpar, logscale=TRUE, cell.labels=cell.label, Lall_mat=lall.mat)



############################################ Sampling Algorithms ####################################################
### Symmetric Proposal distributions

# Proposal_Mu = function(mus.current, pis.current, jumpLength=NULL, K, Smax)
# {
#       mus.next = suppressWarnings(xsample( E = matrix(rep(1,K), nrow=1), F = crossprod(1:Smax, pis.current[-1]), G = rbind(diag(rep(1,K)) , 
#             diag(rep(-1, K))), H = c(rep(0,K),rep(pis.current[1]-1, K)), iter = 2, burnin = 0, type = "mirror",
#             jmp = jumpLength, x0 =mus.current, test = FALSE)$X[2,])
#       mus.next
# }
# Proposal_Mu(mus.current=c(0.5,0.5,0.3), pis.current=c(0.1,0.6,0.2,0.1), jumpLength=0.1, K=3, Smax=3)

# Proposal_Mu = function(mus.current, jumpLength=NULL, K)
# {
#       G = rbind(rep(1,K),diag(rep(1,K)), diag(rep(-1, K)))
#       H = c(1,rep(0,K),rep(-1, K))
#       mus.next = suppressWarnings(xsample(G = G, H = H, 
#             iter = 2, burnin = 0, type = "mirror", jmp = jumpLength, x0 =mus.current, test = FALSE)$X[2,])
#       mus.next
# }

# mus.candidate = Proposal_Mu(mus.current=c(0.5,0.5,0.3), jumpLength=0.1, K=3)
# sum(mus.candidate) 
# pis.candidate = Proposal_Pi(mus.candidate, c(0.1,0.6,0.2,0.1), 0.1, 3, 3)
# print(pis.candidate)   
# crossprod(0:3, pis.candidate)
# Prob.LgivenMuPi(n=inner.iter, K=3, Smax=3, mu=mus.candidate, pis=pis.candidate, sparse.par=GenMatrices(3,3))



# Proposal_Pi = function(mus.current, pis.current, jumpLength=NULL, K, Smax)
# {# pis.current is of length (Smax + 1)
#       E = matrix(c(rep(1,Smax+1), 0:Smax), byrow=TRUE, nrow=2)
#       F = c(1, sum(mus.current))
#       G = rbind(diag(rep(1, Smax+1)), c(-1, rep(0,Smax)))
#       H = c(rep(0, Smax+1), max(mus.current)-1)
#       pis.next = suppressWarnings(xsample( E = E, F = F, G = G, H = H, 
#             iter = 2, burnin = 0, type = "mirror", jmp = jumpLength, x0 =pis.current, test = FALSE)$X[2,])
#       pis.next
# }
# Proposal_Pi(mus.current=c(0.5,0.5,0.3), pis.current=c(0.1,0.6,0.2,0.1), jumpLength=0.1, K=3, Smax=3)

Proposal_MuPi = function(mus.current, pis.current, jumpLength, K, Smax)
{
      E = matrix(c(rep(1,K),-(0:Smax), rep(0,K), rep(1, Smax+1)), byrow=TRUE, nrow=2)
      F = c(0,1)
      G2 = diag(rep(-1, K+1+Smax))
      G2[1:K, K+1] = -1
      G = rbind(diag(rep(1, K+1+Smax)),      #>=0
                G2[1:K,],    #<=1
                c(rep(1,K),rep(0,1+Smax)) )  #sum mu >1 
      H = c(rep(0,K+1+Smax), rep(-1, K),1)
      mupi.next = suppressWarnings(xsample( E = E, F = F, G = G, H = H, 
             iter = 2, burnin = 0, type = "mirror", jmp = jumpLength, x0 =c(mus.current, pis.current), test = FALSE)$X[2,])
      mupi.next
}

# Proposal_MuPi(mus.current=c(0.5,0.5,0.3), pis.current=c(0.1,0.6,0.2,0.1), jumpLength=0.1, K=3, Smax=3)





Proposal_TFpr = function(rate.current, jumpLength=NULL, K)
{
      G = rbind(diag(rep(1,K)), diag(rep(-1, K)))
      H = c(rep(0,K),rep(-1, K))
      rate.next = suppressWarnings(xsample(G = G, H = H, 
            iter = 2, burnin = 0, type = "mirror", jmp = jumpLength, x0 =rate.current, test = FALSE)$X[2,])
      rate.next
}
# Proposal_TFpr(c(0.8,0.8,0.8), 0.05, 3)



# -------------------------------------------------------------------------------------------------------------------------- #
# Random Walk Proposal

RWMH_sampler.v1 = function(K, Smax, MGS, I_GS, MSS, MBS, MBS.ctrl, N.ctrl, 
					   mu.init=NULL, pi.init, ss_tpr.init, bs_tpr.init, bs_fpr.init, DirichPar, BetaPar,
					   iter, inner.iter, burnin, jmp1, jmp2)
{
	ParMatrix = GenMatrices.Sparse(K,Smax)
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      lall.mat = rbind(rep(0,K), ParMatrix$Lmatrix[,1:K])
      cell.labels = apply(lall.mat,1,function(x) paste(as.character(x),collapse=""))

      mus.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      pis.posterior = matrix(NA, nrow=iter+burnin, ncol=Smax+1)
      ss_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_fpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)


      # Initialize
      mus.posterior[1,] = rep(crossprod(1:Smax, pi.init[-1])/K, K)
      pis.posterior[1,] = pi.init
      ss_tpr.posterior[1,] = ss_tpr.init
      bs_tpr.posterior[1,] = bs_tpr.init
      bs_fpr.posterior[1,] = bs_fpr.init

      accept_track = matrix(0, nrow=iter+burnin, ncol=2)
      accept_track[1,] = rep(1, 2)
      alpha_track = accept_track
      
      message("Start block M-H sampling...")

      cell.prob.previous = Prob.LgivenMuPi(n=inner.iter, K=K, Smax=Smax, 
                        mu=mus.posterior[1,], pis=pis.posterior[1,], sparse.par=ParMatrix)

      logP.rate_t1 = Joint.Density.Cpp(K=K, Smax=Smax, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                    latent.cell.probs=cell.prob.previous, 
                    mus=mus.posterior[1,], 
                    pis=pis.posterior[1,], 
                    ss_tpr=ss_tpr.posterior[1,], 
                    bs_tpr=bs_tpr.posterior[1,], 
                    bs_fpr=bs_fpr.posterior[1,], 
                    DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)
      logP.rate_t0 = NA

      for (i in 2:(iter+burnin))
      {
            # record keeping
            print(c(i, accept_track[i-1,]))
            if (i%%25 == 0) 
            {
                  print(round(c(mus.posterior[i-1,]),3))
                  print(round(c(pis.posterior[i-1,]),3))
                  print(round(c(ss_tpr.posterior[i-1,]),3))
                  print("Trailing 200 sample mean:")
                  print(round(apply(pis.posterior[max(1, i-200):(i-1),],2,mean),3))
            }
            if (i%%10 == 0) 
            {
                  print("Average accept rate:")
                  print(apply(accept_track[1:i,],2,mean))
                  print("Trailing 20 accept rate:")
                  print(apply(accept_track[max(1, i-19):i,],2,mean))
            }

            # sample by blocks: (mus, pis), (t/fpr)
            # Mu, Pi -----------------------------------------------------------------------------------------------
            mupi.candidate = Proposal_MuPi(mus.posterior[i-1,], pis.posterior[i-1,], jmp1, K, Smax)
            #print(mupi.candidate)
            mus.candidate = mupi.candidate[1:K]
            pis.candidate = mupi.candidate[-(1:K)]

            cell.prob.candidate = Prob.LgivenMuPi(n=inner.iter, K=K, Smax=Smax, 
                        mu=mus.candidate, pis=pis.candidate, sparse.par=ParMatrix)
            #print(cell.prob.candidate)

            if (is.na(cell.prob.candidate))
            {
                  logP.mupi_t1 = -Inf
            }
            else{
                  logP.mupi_t1 = Joint.Density.Cpp(K=K, Smax=Smax, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.candidate, 
                        mus=mus.candidate, 
                        pis=pis.candidate, 
                        ss_tpr=ss_tpr.posterior[i-1,], bs_tpr=bs_tpr.posterior[i-1,], bs_fpr=bs_fpr.posterior[i-1,], 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)
            }

            logP.mupi_t0 = ifelse(accept_track[i-1,2], logP.rate_t1, logP.rate_t0)
            
            alpha = exp(logP.mupi_t1 - logP.mupi_t0)
            alpha_track[i,1] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  mus.posterior[i,] = mus.candidate
                  pis.posterior[i,] = pis.candidate
                  accept_track[i,1] = 1
                  cell.prob.previous = cell.prob.candidate

            }
            else # reject
            {
                  mus.posterior[i,] = mus.posterior[i-1,]
                  pis.posterior[i,] = pis.posterior[i-1,]
            }

            # (ss_tpr, bs_tpr, bs_fpr) ------------------------------------------------------------------------
            rate.candidate = Proposal_TFpr(c(ss_tpr.posterior[i-1,], bs_tpr.posterior[i-1,], 
                                          bs_fpr.posterior[i-1,]), jmp2, 3*K)

            logP.rate_t0 = ifelse(accept_track[i,1], logP.mupi_t1, logP.mupi_t0)

            logP.rate_t1 = Joint.Density.Cpp(K=K, Smax=Smax, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.previous, 
                        mus=mus.posterior[i,], pis=pis.posterior[i,], 
                        ss_tpr=rate.candidate[1:K], 
                        bs_tpr=rate.candidate[(K+1):(2*K)], 
                        bs_fpr=rate.candidate[(2*K+1):(3*K)], 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)

            alpha = exp(logP.rate_t1 - logP.rate_t0)
            alpha_track[i,2] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  ss_tpr.posterior[i,] = rate.candidate[1:K] 
                  bs_tpr.posterior[i,] = rate.candidate[(K+1):(2*K)]
                  bs_fpr.posterior[i,] = rate.candidate[(2*K+1):(3*K)]
                  
                  accept_track[i,2] = 1
            }
            else # reject
            {
                  ss_tpr.posterior[i,] = ss_tpr.posterior[i-1,] 
                  bs_tpr.posterior[i,] = bs_tpr.posterior[i-1,]
                  bs_fpr.posterior[i,] = bs_fpr.posterior[i-1,]
                  rate.candidate = NA # reset candidate
            }
      }
      posterior = cbind(mus.posterior, pis.posterior, ss_tpr.posterior, bs_tpr.posterior, bs_fpr.posterior)
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:Smax), 
                              paste0("SS_TPR_",1:K), paste0("BS_TPR_",1:K), paste0("BS_FPR_",1:K))

      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}

# betapar = matrix(c(rep(30,6),rep(2,3), c(1.6, 7.5, 3.5, 7.5, 3.5, 1.6, 198, 38, 18)), ncol=2)
# downsample = sample(1:1000,300)
# result = RWMH_sampler.v1(K=3, Smax=3, MGS=M.GS[downsample,], I_GS=I_GS[downsample], MSS=M.SS[downsample,], MBS=M.BS[downsample,], MBS.ctrl=M.BS.ctrl, N.ctrl=1000, 
#                                  mu.init=NULL, pi.init=c(0.1,0.5,0.25,0.15), 
#                                  ss_tpr.init=rep(0.9,3), bs_tpr.init=rep(0.85,3), bs_fpr.init=rep(0.05,3), 
#                                  DirichPar=c(1,4,2,1), BetaPar=betapar,
#                                  iter=400, inner.iter=3000, burnin=50, jmp1=0.01, jmp2=0.005)

# colMeans(result$posterior[-(1:100),])
# stem(result$posterior[-(1:100),3])

RWMH_sampler.v2 = function(K, Smax, HASGS=TRUE, MGS, I_GS, MSS, MBS, MBS.ctrl, N.ctrl, 
                                 mu.init=NULL, pi.init, ss_tpr.init, bs_tpr.init, bs_fpr.init, DirichPar, BetaPar,
                                 iter, inner.iter, burnin, jmp1, jmp2, jmp3, jmp4)
{
      ParMatrix = GenMatrices.Sparse(K,Smax)
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      lall.mat = rbind(rep(0,K), ParMatrix$Lmatrix[,1:K])
      cell.labels = apply(lall.mat,1,function(x) paste(as.character(x),collapse=""))

      mus.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      pis.posterior = matrix(NA, nrow=iter+burnin, ncol=Smax+1)
      ss_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_fpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)


      # Initialize
      mus.posterior[1,] = rep(crossprod(1:Smax, pi.init[-1])/K, K)
      pis.posterior[1,] = pi.init
      ss_tpr.posterior[1,] = ss_tpr.init
      bs_tpr.posterior[1,] = bs_tpr.init
      bs_fpr.posterior[1,] = bs_fpr.init

      accept_track = matrix(0, nrow=iter+burnin, ncol=4)
      accept_track[1,] = rep(1, 4)
      alpha_track = accept_track
      loglik_track = numeric(iter+burnin)
      
      message("Start block M-H sampling...")

      cell.prob.previous = Prob.LgivenMuPi(n=inner.iter, K=K, Smax=Smax, 
                        mu=mus.posterior[1,], pis=pis.posterior[1,], sparse.par=ParMatrix)

      loglik_track[1] = Joint.Density.Cpp(K=K, Smax=Smax, HasGS=HASGS, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                    latent.cell.probs=cell.prob.previous, 
                    mus=mus.posterior[1,], 
                    pis=pis.posterior[1,], 
                    ss_tpr=ss_tpr.posterior[1,], 
                    bs_tpr=bs_tpr.posterior[1,], 
                    bs_fpr=bs_fpr.posterior[1,], 
                    DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)

      for (i in 2:(iter+burnin))
      {
            # record keeping
            print(c(i, accept_track[i-1,],floor(loglik_track[i-1])))
            if (i%%25 == 0) 
            {
                  print(round(c(mus.posterior[i-1,]),3))
                  print(round(c(pis.posterior[i-1,]),3))
                  print(round(c(ss_tpr.posterior[i-1,]),3))
                  print("Trailing 200 sample mean:")
                  print(round(apply(pis.posterior[max(1, i-200):(i-1),],2,mean),3))
            }
            if (i%%10 == 0) 
            {
                  print("Average accept rate:")
                  print(apply(accept_track[1:i,],2,mean))
                  print("Trailing 20 accept rate:")
                  print(apply(accept_track[max(1, i-19):i,],2,mean))
            }

            # sample by blocks: (mus, pis), (t/fpr)
            # Mu, Pi -----------------------------------------------------------------------------------------------
            mupi.candidate = Proposal_MuPi(mus.posterior[i-1,], pis.posterior[i-1,], jmp1, K, Smax)
            #print(mupi.candidate)
            mus.candidate = mupi.candidate[1:K]
            pis.candidate = mupi.candidate[-(1:K)]

            cell.prob.candidate = Prob.LgivenMuPi(n=inner.iter, K=K, Smax=Smax, 
                        mu=mus.candidate, pis=pis.candidate, sparse.par=ParMatrix)
            #print(cell.prob.candidate)

            if (is.na(cell.prob.candidate))
            {
                  logP.mupi_t1 = -Inf
            }
            else{
                  logP.mupi_t1 = Joint.Density.Cpp(K=K, Smax=Smax, HasGS=HASGS, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.candidate, 
                        mus=mus.candidate, 
                        pis=pis.candidate, 
                        ss_tpr=ss_tpr.posterior[i-1,], bs_tpr=bs_tpr.posterior[i-1,], bs_fpr=bs_fpr.posterior[i-1,], 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)
            }

            logP.mupi_t0 = loglik_track[i-1]
            
            alpha = exp(logP.mupi_t1 - logP.mupi_t0)
            alpha_track[i,1] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  mus.posterior[i,] = mus.candidate
                  pis.posterior[i,] = pis.candidate
                  accept_track[i,1] = 1
                  cell.prob.previous = cell.prob.candidate

            }
            else # reject
            {
                  mus.posterior[i,] = mus.posterior[i-1,]
                  pis.posterior[i,] = pis.posterior[i-1,]
            }

            # (ss_tpr) ------------------------------------------------------------------------
            ss_tpr.candidate = Proposal_TFpr(ss_tpr.posterior[i-1,], jmp2, K)

            logP.sstpr_t0 = ifelse(accept_track[i,1], logP.mupi_t1, logP.mupi_t0)

            logP.sstpr_t1 = Joint.Density.Cpp(K=K, Smax=Smax, HasGS=HASGS, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.previous, 
                        mus=mus.posterior[i,], pis=pis.posterior[i,], 
                        ss_tpr=ss_tpr.candidate, 
                        bs_tpr=bs_tpr.posterior[i-1,], 
                        bs_fpr=bs_fpr.posterior[i-1,], 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)

            alpha = exp(logP.sstpr_t1 - logP.sstpr_t0)
            alpha_track[i,2] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  ss_tpr.posterior[i,] = ss_tpr.candidate                   
                  accept_track[i,2] = 1
            }
            else # reject
            {
                  ss_tpr.posterior[i,] = ss_tpr.posterior[i-1,] 
                  ss_tpr.candidate = NA # reset candidate
            }

            # (bs_tpr) ------------------------------------------------------------------------
            bs_tpr.candidate = Proposal_TFpr(bs_tpr.posterior[i-1,], jmp3, K)

            logP.bstpr_t0 = ifelse(accept_track[i,1], logP.sstpr_t1, logP.sstpr_t0)

            logP.bstpr_t1 = Joint.Density.Cpp(K=K, Smax=Smax, HasGS=HASGS, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.previous, 
                        mus=mus.posterior[i,], pis=pis.posterior[i,], 
                        ss_tpr=ss_tpr.posterior[i,], 
                        bs_tpr=bs_tpr.candidate, 
                        bs_fpr=bs_fpr.posterior[i-1,], 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)

            alpha = exp(logP.bstpr_t1 - logP.bstpr_t0)
            alpha_track[i,3] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  bs_tpr.posterior[i,] = bs_tpr.candidate                   
                  accept_track[i,3] = 1
            }
            else # reject
            {
                  bs_tpr.posterior[i,] = bs_tpr.posterior[i-1,] 
                  bs_tpr.candidate = NA # reset candidate
            }


            # (bs_fpr) ------------------------------------------------------------------------
            bs_fpr.candidate = Proposal_TFpr(bs_fpr.posterior[i-1,], jmp4, K)

            logP.bsfpr_t0 = ifelse(accept_track[i,1], logP.bstpr_t1, logP.bstpr_t0)

            logP.bsfpr_t1 = Joint.Density.Cpp(K=K, Smax=Smax, HasGS=HASGS, MGS=MGS, I_GS=I_GS, MSS=MSS, MBS=MBS, Mbs_ctrl=MBS.ctrl, N_ctrl=N.ctrl, 
                        latent.cell.probs=cell.prob.previous, 
                        mus=mus.posterior[i,], pis=pis.posterior[i,], 
                        ss_tpr=ss_tpr.posterior[i,], 
                        bs_tpr=bs_tpr.posterior[i,], 
                        bs_fpr=bs_fpr.candidate, 
                        DirichPar=DirichPar, BetaPars=BetaPar, logscale=TRUE, cell.labels=cell.labels, Lall_mat=lall.mat)

            alpha = exp(logP.bsfpr_t1 - logP.bsfpr_t0)
            alpha_track[i,4] = alpha

            ratio = min(1,alpha); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  bs_fpr.posterior[i,] = bs_fpr.candidate                   
                  accept_track[i,4] = 1
            }
            else # reject
            {
                  bs_fpr.posterior[i,] = bs_fpr.posterior[i-1,] 
                  bs_fpr.candidate = NA # reset candidate
            }

            loglik_track[i] = ifelse(accept_track[i,4], logP.bsfpr_t1, logP.bsfpr_t0)
      }

      posterior = cbind(mus.posterior, pis.posterior, ss_tpr.posterior, bs_tpr.posterior, bs_fpr.posterior)
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:Smax), 
                              paste0("SS_TPR_",1:K), paste0("BS_TPR_",1:K), paste0("BS_FPR_",1:K))

      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),], history.Likelihood = loglik_track[-(1:burnin)])
}
