library(limSolve)
library(doMC)
library(foreach)
library(gtools)
library(Rcpp)
# dependent on ToolBox_PERCH.R

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
# cellLabel = getCellLabel(5,Lmat)

getGSindex = function(MGSi, label)
{
      mgs = paste(as.character(MGSi),collapse="")
      which(label==mgs)
}
# getGSindex(MGS=c(0,1,0), cellLabel)

########################################################################################
################################### Reparameterization #################################
########################################################################################

MuToTheta1 = function(Mu, Theta2, LUmat, MuMat, K, initvalue=NULL){
      Eqt = function(theta1, mu, theta2)
      {
            theta = matrix(c(theta1,theta2),ncol=1)
            # potentials denotes the un-normalized probability
            potentials = exp(LUmat%*%theta)
            # A is the normalizing constant, i.e. sum of all potentials
            A = sum(potentials) + 1
            # k values of mu
            values = as.vector(MuMat%*%potentials) - mu*A
            values
      }

      if (length(initvalue)<K) initvalue = rep(0,K)

      result = nleqslv(initvalue, fn=Eqt, mu = Mu, theta2 = Theta2,
                       method="Broyden",global="dbldog",xscalm="auto", 
                       control = list(maxit=3500,cndtol=1e-15,ftol=1e-10))
      print(result$message)
      return(result)
}

# dmat = AQE.DesignMatrix(3,3)
# MuToTheta1(Mu=c(0.3,0.4,0.5), Theta2=c(-1,-0.5, 0.8), cbind(dmat$Lmat, dmat$Umat), dmat$MuMat, K=3, initvalue=NULL)

########################################################################################
######################## MH within Gibbs with Data Augmentation #########################
########################################################################################

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

GetPosterior = function(K, Smax, HasGS=FALSE, 
                        DataList, # MGS, I_GS, MSS, MBS, MBS.ctrl, N.ctrl, 
                        ParInit, # mu.init=NULL, pi.init, ss_tpr.init, bs_tpr.init, bs_fpr.init, 
                        HyperPar, # DirichPar, BetaPar,
                        Control ) # niter, nburn, jmp1, jmp2, jmp3, jmp4)
{
      # data & parameter input
      ParMatrix = GenMatrices.Sparse(K,Smax)
      LUmat = cbind(ParMatrix$Lmat, ParMatrix$Umat)
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      Lmat.Zero = rbind(rep(0,K), ParMatrix$Lmat)
      cell.labels = apply(Lmat.Zero,1,function(x) paste(as.character(x),collapse=""))

      if (HasGS) {
            I_GS = DataList$I_GS
            MGS = DataList$MGS           
      }
      MSS = DataList$MSS
      MBS = DataList$MBS
      MBS.ctrl = DataList$MBS.ctrl
      X = DataList$X

      N.case = nrow(MBS)
      N.ctrl = nrow(MBS.ctrl)

      iter = Control$iter
      burnin = Control$burnin

      # Initialize
      beta.posterior = array(NA, c(K, ncol(X), iter+burnin))
      theta2.posterior = matrix(NA, nrow=iter+burnin, ncol=choose(K,2))
      ss_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_tpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)
      bs_fpr.posterior = matrix(NA, nrow=iter+burnin, ncol=K)

      beta.posterior[,,1] = ParInit$beta
      theta2.posterior[1,] = ParInit$theta2
      ss_tpr.posterior[1,] = ParInit$ss_tpr
      bs_tpr.posterior[1,] = ParInit$bs_tpr
      bs_fpr.posterior[1,] = ParInit$bs_fpr

      accept_track = matrix(0, nrow=iter+burnin, ncol=4)
      accept_track[1,] = rep(1, 4)
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
