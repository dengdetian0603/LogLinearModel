library(dplyr)
library(magrittr)
library(gtools)
library(doMC)
library(foreach)

#  Simulate Multivariate Bernoulli Distribution
##---------------------------------------------------------------------##
library(MVB)
### QE model K = 4, 10 parameters
n = 5000
theta0 = matrix(c(0.6,0.3,-0.1,-1.2,
                 -0.3,0.2,-3.5,-0.1,-2.0,-10.1),nrow=1)

### design matrix: only has intercept for each parameter
x = matrix(rep(1,1*n),nrow=n)

res = mvb.simu(theta0, x, K = 4, offset = 0)
fitMVB = mvbfit(x, res$response, output = 1);fitMVB

##----------------------------------------------------------------------##
# fit a partially oberved model
L = res$response  # the observed binary vectors
#sum(L[,3])/5000
#hist(apply(L,1,sum))

# function expands observed vector to feature space
expandL = function(l) 
{
      order2 = apply(combn(l,2),2,prod)
      c(l,order2)     
}

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

registerDoMC(4) 
#-- expand the observed vectors
L.expand = foreach(i=1:nrow(L), .combine=rbind) %dopar% {
      expandL(L[i,])
}
#-- get all possible combinations
L.all = as.matrix(AllComb(ncol(L))) 
#-- expand all possible combinations
L.all.exp = foreach(i=1:nrow(L.all), .combine=rbind) %dopar% {
      expandL(L.all[i,])
}

#--------------------- Log likelihood with L -------------------------#
loglik1 = function(theta, L=L.expand, L.All=L.all.exp)
{
      theta = matrix(theta,ncol=1)
      n = nrow(L)
      #if (nrow(theta) != ncol(L)) print("Dimension does not match!")
      part1 = sum(L%*%theta)
      part2 = n*log(sum(exp(L.All%*%theta)))
      part2-part1 # return the negative log likelihood for minimization
}

optim(par=rep(0,10),fn=loglik1,method="BFGS")$par
fitMVB$beta

# Simulate S
mvb.simu(theta0, x, K = 4, offset = 0)$response %>% apply(.,1,sum) %>%
      factor(.,levels=0:4) %>% table(.) %>% as.data.frame(.) -> S.freq


#--------------------- Log likelihood with L and S ---------------------#
loglik2 = function(theta, L=L.expand, L.All=L.all.exp,S=S.freq[,2])
{
      theta = matrix(theta,ncol=1)
      n = nrow(L)
      m = sum(S)
      k = length(S) - 1
      part1 = sum(L%*%theta)
      part2 = (n+m)*log(sum(exp(L.All%*%theta)))
      L.s = apply(L.All[,1:4],1,sum)
      part3 = 0
      for (j in 0:k)
      {
            part3 = part3 + S[j+1]*log(sum(exp(L.All[L.s==j,]%*%theta)))      
      }
      part2-part1-part3 # return the negative log likelihood for minimization
}

optim(par=rep(0,10),fn=loglik2,method="BFGS")$par
# estimates get closer to theta0 but not much.


