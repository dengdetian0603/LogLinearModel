###################################################
# Reparameterization based on meeting on 03/27/2015 
# Mapping from canonical parameters to conditional probability and marginal sum probability
# 
# Updated GammaToTheta.QE() function after meeting on 04/08/2015 
# Updated ThetaToMu(), MuToTheta(), and MuToGamma() functions after meeting on 04/20/2015
###################################################
library(dplyr)
library(gtools)
library(doMC)
library(foreach)
library(MVB)
library(nleqslv)
library(BB)
library(mipfp)

# depend on the AllComb() function
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/LS.model.R")
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

# TODO:
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



