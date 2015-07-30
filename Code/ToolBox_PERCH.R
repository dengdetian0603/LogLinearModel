library(limSolve)
library(doMC)
library(foreach)

# function expands observed vector to feature space
expandL.QE = function(l) 
{
      order2 = apply(combn(l,2),2,prod)
      c(l,order2)     
}
expandL = expandL.QE

# function generates all possible vector of K dimension
AllComb = function(K) 
{
      A = list()
      for (i in 1:K)
      {
            A[[i]]=0:1
      }
      expand.grid(A)
}

# expand to general design matrix of L
expandL.LL = function(l) 
{
      K = length(l)
      L.exp = rep(NA,times=2^K-1)
      L.exp[1:K] = l 
      last_index = K
      for (i in 2:K)
      {
            order.i = apply(combn(l,i),2,prod)
            L.exp[(last_index + 1):(last_index + choose(K,i))] = order.i
            last_index = last_index + choose(K,i)
      }
      L.exp
}


# Calculate Pi based on QE parameter Thetas
QE.ThetaToPi = function(theta, L.all.exp, K)
{
      theta = as.vector(theta)
      pis = c()
      L.s = apply(L.all.exp[,1:K],1,sum)
      Anorm = sum(exp(L.all.exp%*%theta))
      for (s in 0:K)
      {
            pis[s+1] = sum(exp(L.all.exp[L.s==s,]%*%theta))/Anorm
      }
      pis
}

#  General Log-linear model (LL)
## Exchangeable parameterization
PiToTheta = function(Pis)
{
      if (abs(sum(Pis)-1)>1e-7) {print("Please use proper probability.")}
      pi0 = Pis[1]
      pis = Pis[-1]
      K = length(pis)
      Y = matrix(NA,nrow=K,ncol=1)
      for (s in 1:K)
      {
            Y[s,1] = log(pis[s]/pi0/choose(K,s))
      }
      M = matrix(0,nrow=K,ncol=K)
      for (i in 1:K)
      {
            for (j in 1:i)
            {
                  M[i,j] = choose(i,j)
            }
      }
      theta = solve(M,Y)
      theta
}

ThetaToGamma = function(Theta, L.All.expd, k)
{
      theta = matrix(Theta,ncol=1)
      # vector of S_i
      L.s = apply(L.All.expd[,1:k],1,sum)
      # reorder by s
      gamma.order = order(L.s)
      L.All.expd = L.All.expd[gamma.order,]
      L.s = L.s[gamma.order]
      
      
      # potentials denotes the un-normalized probability
      potentials = exp(L.All.expd%*%theta)
      # A is the normalizing constant, i.e. sum of all potentials
      A = sum(potentials)
      
      # k+1 values of pi
      pis = rep(NA,k+1)
      pis[1] = 1/A
      for (j in 1:k)
      {
            pis[j+1] = sum(potentials[L.s==j])/A      
      }
      
      # 2^k -1 values of gamma      
      gamma = rep(NA,2^k-1)
      gamma[2^k-1] = 1
      for (i in 1:(2^k-2))
      {
            s = L.s[i+1]
            gamma[i] = potentials[i+1]/sum(potentials[L.s==s])
      }
      pi.names = paste("pi",0:k,sep="")
      gamma.names = apply(L.All.expd[-1,1:k],1,paste,collapse="")
      
      Gamma = data.frame(Var.name = c(pi.names,gamma.names), values = c(pis,gamma))
      Gamma
}

Equation = function(Theta, Gamma, L.All.expd, k)
{
      fn = rep(NA,2^k-1)
      theta = matrix(Theta,ncol=1)
      L.s = apply(L.All.expd[,1:k],1,sum)
      gamma.order = order(L.s)
      L.All.expd = L.All.expd[gamma.order,]
      L.s = L.s[gamma.order]
      
      potentials = exp(L.All.expd%*%theta)
      A = sum(potentials)
      
      # k unique equations derived from pi
      fn[1] = 1/A - Gamma[1]
      for (j in 2:k)
      {
            fn[j] = sum(potentials[L.s==j-1])/A - Gamma[j]  
      }
      
      # 2^k - 1 - k equations derived from gammas
      j = k + 1; t = 1
      for (s in 1:(k-1))
      {
            for (i in 1:(choose(k,s)-1))
            {
                  fn[j] = potentials[t+i]/sum(potentials[L.s==s]) - Gamma[t+i+k]
                  j = j + 1
            }
            t = t + choose(k,s)
      }
      fn
}

GammaToTheta = function(Gamma0, L.All.expd0, k0, maxNpar0=30)
{
      Equation = function(Theta, Gamma, L.All.expd, k)
      {
            fn = rep(NA,2^k-1)
            theta = matrix(Theta,ncol=1)
            L.s = apply(L.All.expd[,1:k],1,sum)
            gamma.order = order(L.s)
            L.All.expd = L.All.expd[gamma.order,]
            L.s = L.s[gamma.order]
            
            potentials = exp(L.All.expd%*%theta)
            A = sum(potentials)
            As = numeric(k-1)
            
            # k unique equations derived from pi
            fn[1] = 1/A - Gamma[1]
            for (j in 2:k)
            {
                  As[j-1] = sum(potentials[L.s==j-1])
                  fn[j] = As[j-1]/A - Gamma[j]  
            }
            
            # 2^k - 1 - k equations derived from gammas
            j = k + 1; t = 1
            for (s in 1:(k-1))
            {
                  for (i in 1:(choose(k,s)-1))
                  {
                        fn[j] = potentials[t+i]/As[s] - Gamma[t+i+k]
                        j = j + 1
                  }
                  t = t + choose(k,s)
            }
            fn
      } 
      iters = 1
      par0 = rep(0,2^k0-1)
      result = nleqslv(par0, fn=Equation, 
                       Gamma = Gamma0, L.All.expd = L.All.expd0, k = k0,
                       method="Broyden",global="dbldog",xscalm="auto", 
                       control = list(maxit=1500,cndtol=1e-15,ftol=1e-10))
      print(paste("Starting Value Set",iters,":",result$message))
      while(result$termcd>2 & iters<maxNpar0)
      {
            iters = iters + 1
            par0 = rnorm(2^k0-1,sd=0.1)
            result = nleqslv(par0, fn=Equation, 
                             Gamma = Gamma0, L.All.expd = L.All.expd0, k = k0,
                             method="Broyden",global="dbldog",xscalm="auto", 
                             control = list(maxit=1500,cndtol=1e-15,ftol=1e-10))
            print(paste("Starting Value Set",iters,":",result$message))
      } 
      result
}

GammaToTheta.QE = function(Gamma0, L.All.expd0, k0, maxNpar0=30, Full.designMatrix)
{
      if (Full.designMatrix)
      {
            Equation = function(Theta.qe, Gamma, L.All.expd, k)
            {
                  N.qe = k*(k+1)/2
                  fn = rep(NA,N.qe)
                  theta = matrix(c(Theta.qe,rep(0,2^k-1-k*(k+1)/2)),ncol=1)
                  L.s = apply(L.All.expd[,1:k],1,sum)
                  gamma.order = order(L.s)
                  L.All.expd = L.All.expd[gamma.order,]
                  L.s = L.s[gamma.order]
                  
                  potentials = exp(L.All.expd%*%theta)
                  A = sum(potentials)
                  A1 = sum(potentials[L.s==1])
                  A2 = sum(potentials[L.s==2])
                  
                  # 2 unique equations derived from pi
                  fn[1] = 1/A - Gamma[1] # pi0
                  fn[2] = A1/A - Gamma[2] # pi1  
                  
                  # k-1 + (k choose 2)-1 equations derived from gammas
                  j = 3; t = 1
                  for (s in 1:2)
                  {
                        As = ifelse(s==1, A1, A2)
                        for (i in 1:(choose(k,s)-1))
                        {
                              fn[j] = potentials[t+i]/As - Gamma[t+i+k]
                              j = j + 1
                        }
                        t = t + choose(k,s)
                  }
                  fn
            }
      }
      else
      {
            Equation = function(Theta.qe, Gamma, L.All.expd, k)
            {
                  N.qe = k*(k+1)/2
                  fn = rep(NA,N.qe)
                  theta = matrix(Theta.qe,ncol=1)
                  L.s = apply(L.All.expd[,1:k],1,sum)
                  gamma.order = order(L.s)
                  L.All.expd = L.All.expd[gamma.order,]
                  L.s = L.s[gamma.order]
                  
                  potentials = exp(L.All.expd%*%theta)
                  A = sum(potentials)
                  A1 = sum(potentials[L.s==1])
                  A2 = sum(potentials[L.s==2])
                  
                  # 2 unique equations derived from pi
                  fn[1] = 1/A - Gamma[1] # pi0
                  fn[2] = A1/A - Gamma[2] # pi1  
                  
                  # k-1 + (k choose 2)-1 equations derived from gammas
                  j = 3; t = 1
                  for (s in 1:2)
                  {
                        As = ifelse(s==1, A1, A2)
                        for (i in 1:(choose(k,s)-1))
                        {
                              fn[j] = potentials[t+i]/As - Gamma[t+i+k]
                              j = j + 1
                        }
                        t = t + choose(k,s)
                  }
                  fn
            }
      }
      
      iters = 1
      par0 = rep(0,k0*(k0+1)/2)
      result = nleqslv(par0, fn=Equation, 
                       Gamma = Gamma0, L.All.expd = L.All.expd0, k = k0,
                       method="Broyden",global="dbldog",xscalm="auto", 
                       control = list(maxit=3500,cndtol=1e-15,ftol=1e-10))
      print(paste("Starting Value Set",iters,":",result$message))
      while(result$termcd>2 & iters<maxNpar0)
      {
            iters = iters + 1
            par0 = rnorm(k0*(k0+1)/2,sd=0.1)
            result = nleqslv(par0, fn=Equation, 
                             Gamma = Gamma0, L.All.expd = L.All.expd0, k = k0,
                             method="Broyden",global="dbldog",xscalm="auto", 
                             control = list(maxit=3500,cndtol=1e-15,ftol=1e-10))
            print(paste("Starting Value Set",iters,":",result$message))
      } 
      result
}

ThetaToMu = function(Theta, L.All.expd, k)
{
      theta = matrix(Theta,ncol=1)
      # potentials denotes the un-normalized probability
      potentials = exp(L.All.expd%*%theta)
      # A is the normalizing constant, i.e. sum of all potentials
      A = sum(potentials)
      # k values of mu
      mus = rep(NA,k)
      for (j in 1:k)
      {
            mus[j] = sum(potentials[L.All.expd[,j]==1])/A      
      }  
      data.frame(Mu_k = mus)
}

MuToTheta1 = function(Mu, Theta2, L.All.expd, k)
{
      Eqt = function(theta1, mu, theta2, l.All.expd, k0)
      {
            theta = matrix(c(theta1,theta2),ncol=1)
            # potentials denotes the un-normalized probability
            potentials = exp(l.All.expd%*%theta)
            # A is the normalizing constant, i.e. sum of all potentials
            A = sum(potentials)
            # k values of mu
            values = rep(NA,k0)
            for (j in 1:k0)
            {
                  values[j] = mu[j]-sum(potentials[l.All.expd[,j]==1])/A      
            }
            values
      }
      par0 = rep(0,k)
      result = nleqslv(par0, fn=Eqt, mu = Mu, theta2 = Theta2, l.All.expd = L.All.expd ,k0=k,
                       method="Broyden",global="dbldog",xscalm="auto", 
                       control = list(maxit=3500,cndtol=1e-15,ftol=1e-10))
      print(result$message)
      result$x
}

MuTheta2ToGamma = function(Mu, Theta2, L.All.expd, k)
{
      Theta1 = MuToTheta1(Mu, Theta2, L.All.expd, k)
      Gamma = ThetaToGamma(c(Theta1,Theta2), L.All.expd, k)
      Gamma
}


GammaToMu = function(Gamma, L.All.expd, k)
{
      # vector of S_i
      L.s = apply(L.All.expd[,1:k],1,sum)
      # reorder by s
      gamma.order = order(L.s)
      L.All.expd = L.All.expd[gamma.order,]
      L.s = L.s[gamma.order]
      L.All.pos = L.All.expd[-1,]
      
      pis = Gamma[2:(k+1)] # pi 1 to k
      gammas = Gamma[-(1:(k+1))] # gamma 1 to 11...1
      mus = rep(NA,k)
      
      pis.ext = rep(pis, choose(k,1:k))
      prob.joint = gammas*pis.ext
      for (j in 1:k)
      {
            mus[j] = sum(prob.joint[L.All.pos[,j]==1])
      }
      data.frame(Mu_k = mus)
}

# sample measurements from latent state with TPR and FPR
LtoY = function(L,TPR,FPR)
{
      k = ncol(L)
      notL = 1-L
      registerDoMC(4) 
      Y= foreach(i=1:nrow(L), .combine=rbind) %dopar% {
            rbinom(k,1,L[i,]*TPR+notL[i,]*FPR)
      }
      Y
}


# Sample data L based on prior on Pi and Mu, under general log linear model
rDataWithPiMu.LL = function(K, mu.prior, pi.prior, n)
{
      #-- get all possible combinations
      L.all = as.matrix(AllComb(K))[-1,]
      S = apply(L.all,1,sum)
      index = order(S)
      L.all = L.all[index,]
      S = S[index]
      J1 = length(S)
      
      #-- expand all possible combinations
      L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
            expandL.LL(L.all[i,])
      }
      
      PiMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            PiMat[i,] = as.numeric(S==i)
      }
      
      MuMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            MuMat[i,] = as.numeric(L.all.exp[,i]==1)
      }
      
      # number of iter and output should be adaptive according to K
      Phis = xsample( E = rbind(PiMat,MuMat), F = c(pi.prior[-1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                      iter = 6000, output = 1000, type = "cda")
      
      i = sample.int(1000,1)
      print(Phis$X[i,])
      cell.prob = c(1, Phis$X[i,])*pi.prior[1]
      result = rmultinom(1,n,cell.prob)
      label = apply(L.all,1,function(x) paste(as.character(x),collapse=""))
      label = c(paste(rep("0",K),collapse=""),label)
      rownames(result) = label
      result
}


# Sample phi from [phi | mu, pi0]
rPhiGivenPi0Mu.LL = function(K, mu, pi0, burnin=3000, n=1000, method="cda")
{
      #-- get all possible combinations
      L.all = as.matrix(AllComb(K))[-1,]
      S = apply(L.all,1,sum)
      index = order(S)
      L.all = L.all[index,]
      S = S[index]
      J1 = length(S)
      
      #-- expand all possible combinations
      L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
            expandL.LL(L.all[i,])
      }
      
      MuMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            MuMat[i,] = as.numeric(L.all.exp[,i]==1)
      }
      
      # TODO: number of iter and output should be adaptive according to K
      Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi0,mu)/pi0, G = diag(rep(1,J1)), H = rep(0,J1),
                      iter = burnin+n, output = n+1, type = method)
      Phis$X[-1,]
}

# Calculate volume of feasible region by simulation given mu, pi0, and basis of solution space for K = 3
## mu[3] should be the smallest 
GetVol.mupi0 = function(mu, pi0, n)
{
      X0 = rep(0,7)
      X0[4] = (1-pi0-mu[3])/pi0
      X0[5] = (1-pi0-mu[2])/pi0
      X0[6] = (1-pi0-mu[1])/pi0
      X0[7] = (sum(mu)+2*pi0-2)/pi0
      X0 = matrix(X0)
      
      b1 = matrix(c(1,0,0,-1,-1,0,1))
      b2 = matrix(c(0,1,0,-1,0,-1,1))
      b3 = matrix(c(0,0,1,0,-1,-1,1))
      basis = cbind(b1,b2,b3)
      
      A = min((1-mu[2]-pi0)/pi0, (1-mu[3]-pi0)/pi0)
      B = min((1-mu[1]-pi0)/pi0, (1-mu[3]-pi0)/pi0)  
      C = min((1-mu[1]-pi0)/pi0, (1-mu[2]-pi0)/pi0)  
      
      U = runif(n, max=A)
      V = runif(n, max=B)
      W = runif(n, max=C)
      
      UVW = rbind(U,V,W)[,!is.na(W)&!is.na(U)]
      N = ncol(UVW)

      candidate = apply(basis%*%UVW,2,function(x) x + X0)
      index = apply(candidate,2,function(x) sum(x>0&x<(1/pi0-1))>6)
            
      ratio = sum(index)/N
      Vol.candidate = A*B*C
      return(data.frame(Volume = ratio*Vol.candidate, 
                      AcceptRatio = ratio, CandidateSample.n = N, VolSize = Vol.candidate))
}

# convert binary vector y to multinomial accumulative vector
BitoMulti = function(dat, K)
{
      #-- dat should be a matrix wih K columns
      if(!is.matrix(dat))
      {
            dat = matrix(dat, nrow = 1)
      }
      
      #-- get all possible combinations
      L.all = as.matrix(AllComb(K))
      S = apply(L.all,1,sum)
      index = order(S)
      L.all = L.all[index,]
      S = S[index]
      J1 = length(S)
      
      label = apply(L.all,1,function(x) paste(as.character(x),collapse=""))
      dat.lab = apply(dat,1,function(x) paste(as.character(x),collapse=""))
      dat.tab = table(dat.lab)
      
      y = rep(0, 2^K)
      for (i in 1:length(dat.tab))
      {
            j = which(label==names(dat.tab)[i])
            y[j] = dat.tab[i]

      }
      names(y)=label
      y
}

# calculate multinomial density given phi and raw data matrix/vector
dMultinom.phi = function(K, dat, y=NULL , phi, logscale=FALSE)
{
      if (length(y)<1) 
      {
            y = BitoMulti(dat, K)
      }
      A = 1 + sum(phi)
      probs = c(1,phi)/A
      dens = dmultinom(x=y, prob=probs,log=logscale)
      dens
}

# Monte Carlo estimate of P[y|mu,pi0]
dYgivenMuPi0 = function(K, dat, y=NULL, mu, pi0, logscale=FALSE, burnin=3000, n=1000, method="cda")
{
      phis = rPhiGivenPi0Mu.LL(K=K, mu=mu, pi0=pi0, burnin=burnin, n=n, method=method)
      if (length(y)<1)
      {
            y = BitoMulti(dat=dat,K=K)
      }
      X = apply(phis, 1, function(x) dMultinom.phi(K=K, y=y, phi=x, logscale=logscale))
      if (logscale)
      {
            result = log(sum(exp(X))) -log(n)
      }
      else
      {
            result = sum(X)/n
      }
      result
}

# density of the prior p[mu, pi0]
## model 1: [mu] ~ logit.normal(0, Sigma),  [pi0|mu] ~ unif(0, 1-max(mu))
library(mvtnorm)
prior.MuPi0 = function(mu, pi0, Sigma, logscale=FALSE)
{
      K = length(mu)
      beta = log(mu/(1-mu))
      density.mu = dmvnorm(x=beta, mean=rep(0,K), sigma=diag(Sigma), log=logscale)
      density.pi0 = 1/(1-max(mu))
      if (logscale)
      {
            return(density.mu + log(density.pi0))
      }
      else
      {
            return(density.mu*density.pi0)
      }
}

# joint density of p[y, mu, pi0]
density.YMuPi0 = function(K, dat, y=NULL, mu, pi0, SigmaInPrior, logscale=FALSE, inner.burnin, inner.iter, method="cda")
{
      if (length(y)<1)
      {
            y = BitoMulti(dat=dat,K=K)
      }
      f1 = dYgivenMuPi0(K=K, y=y, mu=mu, pi0=pi0, logscale=logscale, burnin=inner.burnin, n=inner.iter, method=method)
      f2 = prior.MuPi0(mu=mu, pi0=pi0, Sigma=SigmaInPrior, logscale=logscale)
      if (logscale)
      {
            return(f1+f2)
      }
      else
      {
            return(f1*f2)
      }
}

## posterior sampler with piror on mu (logitnormal) and pi (unifor(0, 1-max(mu)))
post.mu.pi0 = function(K, mu.init, pi0.init, iter, inner.iter, burnin, inner.burnin, dat, MH.par=c(1,1,1,0.5))
{
      y = BitoMulti(dat=dat,K=K)
      posterior = matrix(NA, nrow=iter, ncol=K+1)
      posterior[1,1:K] = mu.init
      posterior[1,K+1] = pi0.init
      
      beta0.mat = matrix(NA, nrow=iter, ncol=K)
      beta0.mat[1,] = log(mu.init/(1-mu.init))
      
      accept_track = matrix(0, nrow=iter, ncol=K+1)
      accept_track[1,] = rep(1, K+1)
      alpha_track = accept_track
      
      mu.candidate = posterior[1,1:K]
      pi0.candidate = posterior[1,K+1]
      for (i in 2:iter)
      {
            posterior[i,1:K] = mu.candidate
            posterior[i,K+1] = pi0.candidate
            print(c(i,posterior[i,]))
            # sample mu
            for (j in 1:K)
            {
                  # sample logit mu, name it beta0
                  beta0_j.new = beta0.mat[i-1,j] + rnorm(1,0,MH.par[j])
                  mu.candidate[j] = 1/(1+exp(-beta0_j.new))
                  if (mu.candidate[j] >= 1-pi0.candidate | sum(mu.candidate) < 1-pi0.candidate)
                  {
                        mu.candidate[j] = posterior[i,j] # reset candidate
                  }
                  else
                  {
                        log.alpha = density.YMuPi0(K=3, y=y, mu=mu.candidate,pi0=pi0.candidate, Sigma=rep(1.6,3),
                                             logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda") -
                                    density.YMuPi0(K=3, y=y, mu=posterior[i,1:K],pi0=posterior[i,K+1], Sigma=rep(1.6,3),
                                       logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda")
                        alpha_track[i,j] = exp(log.alpha)
                        ratio = min(1,alpha_track[i,j])
                        u = runif(1)
                        if(u <= ratio) # accept
                        {
                              posterior[i,j] = mu.candidate[j]
                              accept_track[i,j] = 1
                        }
                        else # reject
                        {
                              mu.candidate[j] = posterior[i,j] # reset candidate
                        }
                  }
            }
            beta0.mat[i,] = log(posterior[i,1:K]/(1-posterior[i,1:K]))
            # sample pi0
            pi0.candidate = posterior[i-1,K+1]+rnorm(1,0,MH.par[K+1])
            if (pi0.candidate >= 1-max(mu.candidate) | pi0.candidate<=max(0, 1-sum(mu.candidate)))
            {
                  # reject
                  pi0.candidate = posterior[i,K+1]
            }
            else
            {
                  log.alpha = density.YMuPi0(K=3, y=y, mu=mu.candidate,pi0=pi0.candidate, Sigma=rep(1.6,3),
                                             logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda") -
                        density.YMuPi0(K=3, y=y, mu=mu.candidate,pi0=posterior[i,K+1], Sigma=rep(1.6,3),
                                       logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda")
                  alpha_track[i,K+1] = exp(log.alpha)
                  ratio = min(1,alpha_track[i,K+1])
                  u = runif(1)
                  if(u <= ratio) # accept
                  {
                        posterior[i,K+1] = pi0.candidate
                        accept_track[i,K+1] = 1
                  }
                  else # reject
                  {
                        pi0.candidate = posterior[i,K+1]
                  }
            }
      }
      colnames(posterior) = c(paste0("Mu",1:K), "Pi0")
      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}


# ------------------------------------------------------------------------------------------------------------- #
# for mu and pi full vector
GenMatrices = function(K)
{
      L.all = as.matrix(AllComb(K))[-1,]
      S = apply(L.all,1,sum)
      index = order(S)
      L.all = L.all[index,]
      S = S[index]
      J1 = length(S)
      
      #-- expand all possible combinations
      L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
            expandL.LL(L.all[i,])
      }
      
      MuMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            MuMat[i,] = as.numeric(L.all.exp[,i]==1)
      }
      
      PiMat = matrix(NA,nrow=K,ncol=J1)
      for (i in 1:K)
      {
            PiMat[i,] = as.numeric(S==i)
      }
      return(list(MuMat=MuMat, PiMat=PiMat, J1=J1, Lmatrix = L.all.exp))
}

rPhiGivenPiMu.LL = function(K, mu, pi, burnin=30, n=10, method="cda", pars)
{
      # TODO: number of iter and output should be adaptive according to K
      Phis = tryCatch(xsample( E = rbind(pars$PiMat,pars$MuMat), F = c(pi[-1],mu)/pi[1], G = diag(rep(1,pars$J1)), H = rep(0,pars$J1),
                      iter = burnin+n, output = n+1, type = method), error = function(e) "Incompatible constraints.")
      if (is.character(Phis)) return(NA)
      else return( tryCatch(Phis$X[-1,], error=function(e){ print(e); print("pi:");print(pi);print("mu:"); print(mu); return(NA)}))
}

# Monte Carlo estimate of P[y|mu,pi]
density.YgivenMuPi = function(K, dat, y=NULL, mu, pi, logscale=FALSE, burnin=300, n=100, method="cda", parMat)
{
      phis = suppressWarnings(rPhiGivenPiMu.LL(K=K, mu=mu, pi=pi, burnin=burnin, n=n, method=method, pars=parMat))
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
            X = apply(phis, 1, function(x) dMultinom.phi(K=K, y=y, phi=x, logscale=logscale))
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

# density of the prior p[mu, pi]
## model 1: [mu] ~ logit.normal(0, Sigma),  [pi0|mu] ~ unif(0, 1-max(mu)), (pi1,..,piK) ~ Stick Breaking 
prior.MuPi = function(mu, pi, Sigma, Alpha, logscale=FALSE)
{
      K = length(mu)
      beta = log(mu/(1-mu))
      density.mu = dmvnorm(x=beta, mean=rep(0,K), sigma=diag(Sigma), log=logscale)
      density.pi0 = 1/(1-max(mu))
      density.pi1toK = dStickBreak(sticks = pi[-1]/(1-pi[1]), Alpha, logscale)
      if (logscale)
      {
            return(density.mu + log(density.pi0) + density.pi1toK)
      }
      else
      {
            return(density.mu*density.pi0*density.pi1toK)
      }
}

rStickBreak = function(num_weights, alpha) 
{
      betas = rbeta(num_weights, 1, alpha)
      remaining_stick_lengths = c(1, cumprod(1 - betas))[1:num_weights]
      weights = remaining_stick_lengths * betas
      weights
}

dStickBreak = function(sticks, alpha, logscale=FALSE)
{
      k = length(sticks)
      w = c(1,1-cumsum(sticks))[1:k]
      betas = sticks/w
      if (logscale) return(sum(dbeta(betas[-k], 1, alpha, log=logscale)))
      else return(prod(dbeta(betas[-k], 1, alpha, log=logscale)))
}


# joint density of p[y, mu, pi]
density.YMuPi = function(K, dat, y=NULL, mu, pi, SigmaInPrior, AlphaInPrior, logscale=FALSE, inner.burnin, inner.iter, method="cda", ParMat)
{
      if (length(y)<1)
      {
            y = BitoMulti(dat=dat,K=K)
      }
      f1 = density.YgivenMuPi(K=K, y=y, mu=mu, pi=pi, logscale=logscale, burnin=inner.burnin, n=inner.iter, method=method, parMat=ParMat)
      f2 = prior.MuPi(mu=mu, pi=pi, Sigma=SigmaInPrior, Alpha=AlphaInPrior, logscale=logscale)
      if (logscale)
      {
            return(f1+f2)
      }
      else
      {
            return(f1*f2)
      }
}

## posterior sampler with piror on mu (logitnormal) and pi (unifor(0, 1-max(mu)))
# TODO!!!
post.mu.pi.ByBlock = function(K, mu.init=NULL, pi.init, iter, inner.iter, burnin, inner.burnin, dat, MH.par=c(1,1,1,0.5), ParMatrix, prior.alpha)
{
      y = BitoMulti(dat=dat,K=K)
      posterior = matrix(NA, nrow=iter, ncol=2*K+1)
      posterior[1,(K+1):(2*K+1)] = pi.init
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      
      tmp = xsample( E = PiMat, F = pi.init[-1]/pi.init[1], G = diag(rep(1,J1)), H = rep(0,J1),
                     iter = 10, output = 5, type = "cda")
      mu.sample = MuMat%*%t(tmp$X)*pi.init[1]  
      posterior[1,1:K] = mu.sample[,2] # TODO: random init?
      
      # beta0.mat = matrix(NA, nrow=iter, ncol=K)
      # beta0.mat[1,] = log(posterior[1,1:K]/(1-posterior[1,1:K]))
      
      accept_track = matrix(0, nrow=iter, ncol=2)
      accept_track[1,] = rep(1, 2)
      alpha_track = accept_track
      
      mu.candidate = posterior[1,1:K]
      pi.candidate = posterior[1,(K+1):(2*K+1)]
      message("Start block M-H sampling...")
      for (i in 2:iter)
      {
            posterior[i,1:K] = mu.candidate
            posterior[i,(K+1):(2*K+1)] = pi.candidate
            print(i)
            if (i%%25 == 0) print(c(i,posterior[i,]))
            # sample mu
            
            tmp = xsample( E = PiMat, F = pi.candidate[-1]/pi.candidate[1], G = diag(rep(1,J1)), H = rep(0,J1),
                           iter = 50, output = 40, type = "cda")
            mu.sample = MuMat%*%t(tmp$X)*pi.candidate[1]  
            mu.candidate = mu.sample[,sample(1:40,1)]
            
            log.alpha = density.YMuPi(K=K, y=y, mu=mu.candidate,pi=pi.candidate, SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix) -
                  density.YMuPi(K=K, y=y, mu=posterior[i,1:K],pi=pi.candidate, SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix)
            alpha_track[i,1] = exp(log.alpha)
            ratio = min(1,alpha_track[i,1])
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
#print("checkmark1")
            pi0.candidate = posterior[i-1,K+1]+rnorm(1,0,MH.par[K+1])
            if (pi0.candidate >= 1-max(mu.candidate) | pi0.candidate<=max(0, 1-sum(mu.candidate)))
            {
                  # reject
                  pi.candidate = posterior[i,(K+1):(2*K+1)]
            }
            else
            {
#print("checkmark2")
                  Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi0.candidate,posterior[i,1:K])/pi0.candidate, G = diag(rep(1,J1)), H = rep(0,J1),
                                  iter = 50, output = 40, type = "cda")
                  pi.sample = PiMat%*%t(Phis$X)*pi0.candidate
                  pi.candidate = c(pi0.candidate, pi.sample[,sample(1:40,1)])
                  
                  log.alpha = density.YMuPi(K=K, y=y, mu=mu.candidate,pi=pi.candidate, SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                            logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix) -
                        density.YMuPi(K=K, y=y, mu=mu.candidate,pi=posterior[i,(K+1):(2*K+1)], SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                      logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix)
                  alpha_track[i,2] = exp(log.alpha)
                  ratio = min(1,alpha_track[i,2])
                  u = runif(1)
#print("checkmark3")
                  if(u <= ratio) # accept
                  {
                        posterior[i, (K+1):(2*K+1)] = pi.candidate
                        accept_track[i,2] = 1
                  }
                  else # reject
                  {
                        pi.candidate = posterior[i,(K+1):(2*K+1)]
                  }
            }
#print("checkmark4")
      }
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:K))
      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}


post.mu.pi = function(K, mu.init=NULL, pi.init, iter, inner.iter, burnin, inner.burnin, dat, MH.par=c(1,1,1,0.5), ParMatrix, prior.alpha)
{
      y = BitoMulti(dat=dat,K=K)
      posterior = matrix(NA, nrow=iter, ncol=2*K+1)
      posterior[1,(K+1):(2*K+1)] = pi.init
      PiMat = ParMatrix$PiMat
      MuMat = ParMatrix$MuMat
      J1 = ParMatrix$J1
      
      tmp = xsample( E = PiMat, F = pi.init[-1]/pi.init[1], G = diag(rep(1,J1)), H = rep(0,J1),
                      iter = 50, output = 10, type = "cda")
      mu.sample = MuMat%*%t(tmp$X)*pi.init[1]  
      posterior[1,1:K] = mu.sample[,2] # TODO: random init?
      
      beta0.mat = matrix(NA, nrow=iter, ncol=K)
      beta0.mat[1,] = log(posterior[1,1:K]/(1-posterior[1,1:K]))
      
      accept_track = matrix(0, nrow=iter, ncol=K+1)
      accept_track[1,] = rep(1, K+1)
      alpha_track = accept_track
      
      mu.candidate = posterior[1,1:K]
      pi.candidate = posterior[1,(K+1):(2*K+1)]
      for (i in 2:iter)
      {
            posterior[i,1:K] = mu.candidate
            posterior[i,(K+1):(2*K+1)] = pi.candidate
            
            print(c(i,posterior[i,]))
            # sample mu
            for (j in 1:K)
            {
                  # sample logit mu, name it beta0
                  beta0_j.new = beta0.mat[i-1,j] + rnorm(1,0,MH.par[j])
                  mu.candidate[j] = 1/(1+exp(-beta0_j.new))
                  if (mu.candidate[j] >= 1-pi.candidate[1] | sum(mu.candidate) < 1-pi.candidate[1])
                  {
                        mu.candidate[j] = posterior[i,j] # reset candidate
                  }
                  else
                  {
                        log.alpha = density.YMuPi(K=K, y=y, mu=mu.candidate,pi=pi.candidate, SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                                   logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix) -
                              density.YMuPi(K=K, y=y, mu=posterior[i,1:K],pi=posterior[i,(K+1):(2*K+1)], SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                             logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix)
                        alpha_track[i,j] = exp(log.alpha)
                        ratio = min(1,alpha_track[i,j])
                        u = runif(1)
                        if(u <= ratio) # accept
                        {
                              posterior[i,j] = mu.candidate[j]
                              accept_track[i,j] = 1
                        }
                        else # reject
                        {
                              mu.candidate[j] = posterior[i,j] # reset candidate
                        }
                  }
            }
            beta0.mat[i,] = log(posterior[i,1:K]/(1-posterior[i,1:K]))
            # sample pi0
            pi0.candidate = posterior[i-1,K+1]+rnorm(1,0,MH.par[K+1])
            if (pi0.candidate >= 1-max(mu.candidate) | pi0.candidate<=max(0, 1-sum(mu.candidate)))
            {
                  # reject
                  pi.candidate = posterior[i,(K+1):(2*K+1)]
            }
            else
            {
                  Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi0.candidate,posterior[i,1:K])/pi0.candidate, G = diag(rep(1,J1)), H = rep(0,J1),
                                  iter = 50, output = 40, type = "cda")
                  pi.sample = PiMat%*%t(Phis$X)*pi0.candidate
                  pi.candidate = c(pi0.candidate, pi.sample[,sample(1:40,1)])
                  
                  log.alpha = density.YMuPi(K=K, y=y, mu=mu.candidate,pi=pi.candidate, SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                             logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix) -
                        density.YMuPi(K=K, y=y, mu=mu.candidate,pi=posterior[i,(K+1):(2*K+1)], SigmaInPrior=rep(1.6,K), AlphaInPrior=prior.alpha,
                                       logscale=TRUE, inner.burnin=inner.burnin, inner.iter=inner.iter, method="cda", ParMat=ParMatrix)
                  alpha_track[i,K+1] = exp(log.alpha)
                  ratio = min(1,alpha_track[i,K+1])
                  u = runif(1)
                  if(u <= ratio) # accept
                  {
                        posterior[i, (K+1):(2*K+1)] = pi.candidate
                        accept_track[i,K+1] = 1
                  }
                  else # reject
                  {
                        pi.candidate = posterior[i,(K+1):(2*K+1)]
                  }
            }
      }
      colnames(posterior) = c(paste0("Mu_",1:K), paste0("Pi_",0:K))
      list(posterior = posterior[-(1:burnin),], history.alpha = alpha_track[-(1:burnin),], 
           history.accept = accept_track[-(1:burnin),])
}



