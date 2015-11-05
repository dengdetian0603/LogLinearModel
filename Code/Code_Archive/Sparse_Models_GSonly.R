library(limSolve)
library(doMC)
library(foreach)
library(gtools)
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
      return(list(MuMat=MuMat, PiMat=PiMat, J1=J1))
}


rPhiGivenPiMu.Sparse = function(K, Smax, mu, pis, burnin=30, n=10, method="cda", pars)
{
      # TODO: number of iter and output should be adaptive according to K
      pis = pis[1:(Smax+1)]
      if (abs(sum(pis)-1)>1e-9)
      {
            print(sum(pis))
            message("pi vector is not Sparse.")
            return(NA)
      }
      Phis = tryCatch(xsample( E = rbind(pars$PiMat,pars$MuMat), F = c(pis[-1],mu)/pis[1], G = diag(rep(1,pars$J1)), H = rep(0,pars$J1),
                      iter = n, burnin = burnin, type = method, test=FALSE), error = function(e) "Incompatible constraints.")
      if (is.character(Phis)) return(NA)
      else return( tryCatch(Phis$X[-1,], error=function(e){ print(e); print("pi:");print(pis);print("mu:"); print(mu); return(NA)}))
}

#rPhiGivenPiMu.Sparse(K=5,Smax=3,mu=TrueMu.sparse, pis=TruePi.sparse, method="mirror",pars=Lmat)


dMultinom.phi.Sparse = function(K, dat, y=NULL , phi, logscale=FALSE)
{
      if (length(y)<1) 
      {
            y = BitoMulti(dat, K)
      }
      J = length(phi)
      if (J < 2^K - 1)
      {
            phi[(J+1):(2^K-1)]=0
      }
      A = 1 + sum(phi)
      probs = c(1,phi)/A
      dens = dmultinom(x=y, prob=probs,log=logscale)
      dens
}


density.YgivenMuPi.Sparse = function(K, Smax, dat, y=NULL, mu, pis, logscale=FALSE, burnin=300, n=100, method="mirror", parMat)
{
      phis = suppressWarnings(rPhiGivenPiMu.Sparse(K=K, Smax=Smax, mu=mu, pis=pis, burnin=burnin, n=n, method=method, pars=parMat))
      if (length(phis)<=1) 
      {
            warning("Incompatible Mu,Pi value, density value is considered 0.")
            if (logscale)
            {return(log(0))}
            else
            {return(0)}
      }
      else
      {
            if (length(y)<1)
            {
                  y = BitoMulti(dat=dat,K=K)
            }
            J = ncol(phis)
            R = nrow(phis)
            phis = cbind(phis, matrix(0, nrow=R, ncol=2^K-1-J))
            X = apply(phis, 1, function(x) dMultinom.phi.Sparse(K=K, y=y, phi=x, logscale=logscale))
            if (logscale)
            {
                  result = log(sum(exp(X))) -log(n)
            }
            else
            {
                  result = sum(X)/n
            }
            return(result)
      }
}

#density.YgivenMuPi.Sparse(5,3,dat=dat, mu=TrueMu.sparse, pis=TruePi.sparse, logscale=TRUE, parMat=Lmat_S)

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

density.YMuPi.Sparse = function(K, Smax, dat, y=NULL, mu, pis, AlphaInPrior, logscale=FALSE, inner.burnin, inner.iter, method="cda", ParMat)
{
      if (length(y)<1)
      {
            y = BitoMulti(dat=dat,K=K)
      }
      f1 = density.YgivenMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu, pis=pis, logscale=logscale, burnin=inner.burnin, n=inner.iter, method=method, parMat=ParMat)
      f2 = prior.MuPi.Sparse(mu=mu, pis=pis, Alpha=AlphaInPrior, Smax=Smax, logscale=logscale)
      if (logscale)
      {
            return(f1+f2)
      }
      else
      {
            return(f1*f2)
      }
}

#prior.MuPi.Sparse(mu=TrueMu.sparse, pis=TruePi.sparse, Alpha=c(1,4,2,1), Smax=3, logscale=TRUE)
#density.YMuPi.Sparse(K=5, Smax=3, dat=dat[1:100,], mu=TrueMu.sparse, pis=TruePi.sparse, AlphaInPrior=c(1,4,2,1), logscale=TRUE, 
#    inner.burnin=50, inner.iter=2000, method="mirror", ParMat=Lmat_S)

# -------------------------------------------------------------------------------------------------------------------------- #
# Random Walk Proposal
post.mu.pi.Sparse.v1 = function(K, Smax, mu.init=NULL, pi.init, iter, inner.iter, burnin, inner.burnin, dat, 
      MH.sigmaOfpi0=0.1, MH.sigmaOfpi.rest=rep(0.3, Smax-2), MH.sigmaOfmu=rep(1,K-1), ParMatrix, prior.alpha, 
      densityMethod = "mirror", Ncore=detectCores())
{
      y = BitoMulti(dat=dat,K=K)
      posterior = matrix(NA, nrow=iter+burnin, ncol=K+Smax+1)
      posterior[1,(K+1):(K+Smax+1)] = pi.init
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      
      # Initialize
      tmp = xsample( E = matrix(rep(1,K), nrow=1), F = crossprod(1:Smax, pi.init[-1]), G = rbind(diag(rep(1,K)) , 
            diag(rep(-1, K))), H = c(rep(0,K),rep(pi.init[1]-1, K)), iter = 100, burnin = 50, type = "mirror")$X
      mu.sample = tmp[sample(1:100,1),]
      posterior[1,1:K] = mu.sample 

      # beta0.mat = matrix(NA, nrow=iter, ncol=K)
      # beta0.mat[1,] = log(posterior[1,1:K]/(1-posterior[1,1:K]))
      
      accept_track = matrix(0, nrow=iter+burnin, ncol=2)
      accept_track[1,] = rep(1, 2)
      alpha_track = accept_track
      
      mu.candidate = posterior[1,1:K]
      pi.candidate = posterior[1,(K+1):(Smax+K+1)]
      message("Start block M-H sampling...")

      for (i in 2:(iter+burnin))
      {
            posterior[i,1:K] = mu.candidate
            posterior[i,(K+1):(Smax+K+1)] = pi.candidate
            print(c(i, accept_track[i-1,]))
            if (i%%25 == 0) 
            {
                  print(round(posterior[i,],3))
                  print("Trailing 200 sample mean:")
                  print(round(apply(posterior[max(1, i-199):i,],2,mean),3))
            }
            if (i%%10 == 0) 
            {
                  print("Average accept rate:")
                  print(apply(accept_track[1:i,],2,mean))
                  print("Trailing 20 accept rate:")
                  print(apply(accept_track[max(1, i-19):i,],2,mean))
            }

            # sample mu by block
            mu.candidate[1:(K-1)] = inv.logit(logit(posterior[i-1,1:(K-1)]) + rnorm(K-1, 0, MH.sigmaOfmu))
            mu.candidate[K] = crossprod(1:Smax, posterior[i-1,(K+2):(Smax+K+1)]) - sum(mu.candidate[1:(K-1)])
            
            if (mu.candidate[K]>0 & mu.candidate[K]<1)
            {
                  logJointDensity = foreach(lik = 1:Ncore, .combine=c) %dopar% {
                        if (lik < Ncore/2 + 1){
                              density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=pi.candidate, AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter/2, method=densityMethod, ParMat=ParMatrix)
                              } else {
                              density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=posterior[i,1:K],pis=pi.candidate, AlphaInPrior=prior.alpha,
                                logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter/2, method=densityMethod, ParMat=ParMatrix)      
                              }
                  }
                  #print(logJointDensity)
                  log.alpha = mean(logJointDensity[1:(Ncore/2)] - logJointDensity[(Ncore/2 + 1):Ncore]) + 
                  sum(log(deriv.logit(posterior[i,1:(K-1)]))) - sum(log(deriv.logit(mu.candidate[1:(K-1)])))
            } else {
                  log.alpha = -Inf
            }
            alpha_track[i,1] = exp(log.alpha)
            ratio = min(1,alpha_track[i,1]); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  posterior[i,1:K] = mu.candidate
                  accept_track[i,1] = 1
            }
            else # reject
            {
                  mu.candidate = posterior[i,1:K] # reset candidate
            }

            #beta0.mat[i,] = log(posterior[i,1:K]/(1-posterior[i,1:K]))
            # sample pi0
            pi0.candidate = inv.logit(logit(posterior[i-1,K+1])+rnorm(1,0,MH.sigmaOfpi0))
            if (pi0.candidate >= 1-max(mu.candidate) | pi0.candidate<=max(0, 1-sum(mu.candidate)))
            {
                  # reject
                  pi.candidate = posterior[i,(K+1):(Smax+K+1)]
            }
            else
            {
                  pi1toSmax_2 = inv.logit(logit(posterior[i-1,(K+2):(K+Smax-1)]) + rnorm(Smax-2, 0, MH.sigmaOfpi.rest))
                  b=numeric()
                  b[1] = 1 - pi0.candidate - sum(pi1toSmax_2)
                  b[2] = sum(mu.candidate) - crossprod(1:(Smax-2), pi1toSmax_2)
                  matA = matrix(c(1, Smax-1, 1, Smax), nrow=2)
                  pi.rest = solve(matA,b)

                  pi.candidate = c(pi0.candidate, pi1toSmax_2, pi.rest)

                  if (prod(pi.rest>0 & pi.rest<1)>0)
                  {
                        logJointDensity2 = foreach(lik = 1:Ncore, .combine=c) %dopar% {
                              if (lik < Ncore/2 + 1){
                                    density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=pi.candidate, AlphaInPrior=prior.alpha,
                                            logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method=densityMethod, ParMat=ParMatrix)
                                    } else {
                                    density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=posterior[i,(K+1):(Smax+K+1)], AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method=densityMethod, ParMat=ParMatrix)      
                                    }
                        }
                        #print(logJointDensity2)
                        #print(mean(logJointDensity2[1:(Ncore/2)] - logJointDensity2[(Ncore/2 + 1):Ncore]))
                        #print(sum(log(deriv.logit(posterior[i,(K+1):(Smax+K+1)]))))
                        #print(sum(log(deriv.logit(pi.candidate))))
                        #print(pi.candidate)
                        log.alpha = mean(logJointDensity2[1:(Ncore/2)] - logJointDensity2[(Ncore/2 + 1):Ncore]) + 
                        sum(log(deriv.logit(posterior[i,(K+1):(Smax+K-1)]))) - sum(log(deriv.logit(pi.candidate[1:(Smax-1)])))
                  } else {
                        log.alpha = -Inf
                  }
                  alpha_track[i,2] = exp(log.alpha)
                  ratio = min(1,alpha_track[i,2]); #print(ratio)
                  u = runif(1)

                  if(u <= ratio) # accept
                  {
                        posterior[i, (K+1):(Smax+K+1)] = pi.candidate
                        accept_track[i,2] = 1
                  }
                  else # reject
                  {
                        pi.candidate = posterior[i,(K+1):(Smax+K+1)]
                  }
            }

      }
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:Smax))
      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}

# Uniform Proposal
post.mu.pi.Sparse.v2 = function(K, Smax, mu.init=NULL, pi.init, iter, inner.iter, burnin, inner.burnin, dat, 
      MH.sigmaOfpi0=0.1, MH.sigmaOfpi.rest=rep(0.3, Smax-2), MH.sigmaOfmu=rep(1,K-1), ParMatrix, prior.alpha, 
      densityMethod = "mirror", Ncore=detectCores())
{
      y = BitoMulti(dat=dat,K=K)
      posterior = matrix(NA, nrow=iter+burnin, ncol=K+Smax+1)
      posterior[1,(K+1):(K+Smax+1)] = pi.init
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      
      # Initialize
      tmp = xsample( E = matrix(rep(1,K), nrow=1), F = crossprod(1:Smax, pi.init[-1]), G = rbind(diag(rep(1,K)) , 
            diag(rep(-1, K))), H = c(rep(0,K),rep(pi.init[1]-1, K)), iter = 100, burnin = 50, type = "mirror")$X
      mu.sample = tmp[sample(1:100,1),]
      posterior[1,1:K] = mu.sample 

      # beta0.mat = matrix(NA, nrow=iter, ncol=K)
      # beta0.mat[1,] = log(posterior[1,1:K]/(1-posterior[1,1:K]))
      
      accept_track = matrix(0, nrow=iter+burnin, ncol=2)
      accept_track[1,] = rep(1, 2)
      alpha_track = accept_track
      
      mu.candidate = posterior[1,1:K]
      pi.candidate = posterior[1,(K+1):(Smax+K+1)]
      message("Start block M-H sampling...")

      for (i in 2:(iter+burnin))
      {
            posterior[i,1:K] = mu.candidate
            posterior[i,(K+1):(Smax+K+1)] = pi.candidate
            print(c(i, accept_track[i-1,]))
            if (i%%25 == 0) 
            {
                  print(round(posterior[i,],3))
                  print("Trailing 200 sample mean:")
                  print(round(apply(posterior[max(1, i-199):i,],2,mean),3))
            }
            if (i%%10 == 0) 
            {
                  print("Average accept rate:")
                  print(apply(accept_track[1:i,],2,mean))
                  print("Trailing 20 accept rate:")
                  print(apply(accept_track[max(1, i-19):i,],2,mean))
            }

            # sample mu by block            
            if (accept_track[i-1,2]>0.5 | mean(accept_track[max(1, i-5):(i-1),1]) < 0.05)
            {
            	mu.tmp = xsample( E = matrix(rep(1,K), nrow=1), F = crossprod(1:Smax, posterior[i-1,(K+2):(Smax+K+1)]), 
            				G = rbind(diag(rep(1,K)), diag(rep(-1, K))), 
            				H = c(rep(0,K),rep(posterior[i-1,(K+1)]-1, K)), 
            				iter = 100, burnin = 50, type = "mirror")$X
            }
      	mu.candidate = mu.tmp[sample(2:100,1),]	

		logJointDensity = foreach(lik = 1:Ncore, .combine=c) %dopar% {
                        if (lik < Ncore/2 + 1){
                              density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=pi.candidate, AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter/2, method=densityMethod, ParMat=ParMatrix)
                        } else {
                              density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=posterior[i,1:K],pis=pi.candidate, AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter/2, method=densityMethod, ParMat=ParMatrix)      
                        }
                  }
            #print(logJointDensity)
            log.alpha = mean(logJointDensity[1:(Ncore/2)] - logJointDensity[(Ncore/2 + 1):Ncore])
 
            alpha_track[i,1] = exp(log.alpha)
            ratio = min(1,alpha_track[i,1]); #print(ratio)
            u = runif(1)
            if(u <= ratio) # accept
            {
                  posterior[i,1:K] = mu.candidate
                  accept_track[i,1] = 1
            }
            else # reject
            {
                  mu.candidate = posterior[i,1:K] # reset candidate
            }

            # sample pi0
            pi.tmp = rdirichlet(1, alpha=prior.alpha)
            pi0.candidate = pi.tmp[1] #runif(1, min = max(0, 1-sum(mu.candidate)), max = 1-max(mu.candidate))
            pi1toSmax_2 = pi.tmp[2:(Smax-1)] # runif(Smax-2, min = 0, max = 1-pi0.candidate)
            		  		         #inv.logit(logit(posterior[i-1,(K+2):(K+Smax-1)]) + rnorm(Smax-2, 0, MH.sigmaOfpi.rest))
            b=numeric()
            b[1] = 1 - pi0.candidate - sum(pi1toSmax_2)
            b[2] = sum(mu.candidate) - crossprod(1:(Smax-2), pi1toSmax_2)
            matA = matrix(c(1, Smax-1, 1, Smax), nrow=2)
            pi.rest = solve(matA,b)

            pi.candidate = c(pi0.candidate, pi1toSmax_2, pi.rest)

            if (prod(pi.rest>0 & pi.rest<1)>0)
            {
            	logJointDensity2 = foreach(lik = 1:Ncore, .combine=c) %dopar% {
            		if (lik < Ncore/2 + 1){
            			density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=pi.candidate, AlphaInPrior=prior.alpha,
            				logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method=densityMethod, ParMat=ParMatrix)
            		} else {
            			density.YMuPi.Sparse(K=K, Smax=Smax, y=y, mu=mu.candidate,pis=posterior[i,(K+1):(Smax+K+1)], AlphaInPrior=prior.alpha,
            				logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method=densityMethod, ParMat=ParMatrix)      
            		}
            	}
                  #print(logJointDensity2)
                  #print(mean(logJointDensity2[1:(Ncore/2)] - logJointDensity2[(Ncore/2 + 1):Ncore]))
                  #print(sum(log(deriv.logit(posterior[i,(K+1):(Smax+K+1)]))))
                  #print(sum(log(deriv.logit(pi.candidate))))
                  #print(pi.candidate)
                  log.alpha = mean(logJointDensity2[1:(Ncore/2)] - logJointDensity2[(Ncore/2 + 1):Ncore])
                  
            } else {
                  log.alpha = -Inf
            }
            alpha_track[i,2] = exp(log.alpha)
            ratio = min(1,alpha_track[i,2]); #print(ratio)
            u = runif(1)

            if(u <= ratio) # accept
            {
            	posterior[i, (K+1):(Smax+K+1)] = pi.candidate
                 	accept_track[i,2] = 1
            } else { # reject
                  pi.candidate = posterior[i,(K+1):(Smax+K+1)]
            }

      }
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:Smax))
      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}