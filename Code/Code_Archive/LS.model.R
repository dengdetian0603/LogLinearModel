###################################################
# Reparameterization based on meeting on Mar 9th 
# L be fulling observed binary vector, n observations.
# S be the marginal sums of L that are not fully observed, m observations.
# Combine L and S, solve MLE, to mimick the behavior of model with prior on S
###################################################

library(dplyr)
library(magrittr)
library(gtools)
library(doMC)
library(foreach)
library(MVB)

#  Simulate Multivariate Bernoulli Distribution
##---------------------------------------------------------------------##
### QE model K = 4, 10 parameters
#n = 500
#theta0 = matrix(c(0.6,0.3,-0.1,-1.2,
#                 -0.3,0.2,-3.5,-0.1,-2.0,-10.1),nrow=1)

### design matrix: only has intercept for each parameter
#x = matrix(rep(1,1*n),nrow=n)

#res = mvb.simu(theta0, x, K = 4, offset = 0)
#fitMVB = mvbfit(x, res$response, output = 1);fitMVB

##----------------------------------------------------------------------##

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

mvb.exch.sim = function(pis,n)
{
      k = length(pis) - 1
      s = rmultinom(1,n,pis)
      pool = AllComb(k)
      pool.s = apply(pool,1,sum)
      data = matrix(NA,nrow=n,ncol=k)
      current.row = 0
      for (j in 1:(k+1))
      {
            if(s[j,1]>0)
            {
                  index = sample(1:choose(k,j-1),s[j,1],replace=TRUE)
                  subclass = pool[pool.s==j-1,]
                  tmp = as.matrix(subclass[index,1:k])
                  data[(current.row+1):(current.row+s[j,1]),1:k] = tmp
                  current.row = current.row + s[j,1]
            }
      }
      data
}

#mvb.exch.sim(c(0.3,0.4,0.2,0.05,0.05),30)

#registerDoMC(4) 
#-- expand the observed vectors
#L.expand = foreach(i=1:nrow(L), .combine=rbind) %dopar% {
#      expandL(L[i,])
#}
#-- get all possible combinations
#L.all = as.matrix(AllComb(ncol(L))) 
#-- expand all possible combinations
#L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
#      expandL(L.all[i,])
#}

#--------------------- Log likelihood with L -------------------------#
loglik1 = function(theta, L, L.All)
{
      theta = matrix(theta,ncol=1)
      n = nrow(L)
      #if (nrow(theta) != ncol(L)) print("Dimension does not match!")
      part1 = sum(L%*%theta)
      part2 = n*log(sum(exp(L.All%*%theta)))
      part2-part1 # return the negative log likelihood for minimization
}



# Simulate S
#mvb.simu(theta0, x, K = 4, offset = 0)$response %>% apply(.,1,sum) %>%
#      factor(.,levels=0:4) %>% table(.) %>% as.data.frame(.) -> S.freq


#--------------------- Log likelihood with L and S ---------------------#
loglik2 = function(theta, L, L.All,S)
{
      theta = matrix(theta,ncol=1)
      n = nrow(L)
      m = sum(S)
      k = length(S) - 1
      part1 = sum(L%*%theta)
      part2 = (n+m)*log(sum(exp(L.All%*%theta)))
      L.s = apply(L.All[,1:k],1,sum)
      part3 = 0
      for (j in 0:k)
      {
            part3 = part3 + S[j+1]*log(sum(exp(L.All[L.s==j,]%*%theta)))      
      }
      part2-part1-part3 # return the negative log likelihood for minimization
}

#optim(par=rep(0,10),fn=loglik2,method="BFGS")$par

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
