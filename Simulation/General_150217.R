library(dplyr)
library(magrittr)
library(gtools)
library(doMC)
library(foreach)

##--------------------------------------------------------------------------##
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

PiToTheta(c(0.3,0.4,0.2,0.1))
rdirichlet(1,c(3,4,2,0.5,0.5)) %>% PiToTheta(.)

# simulate
registerDoMC(2)
theta.prior = foreach(i=1:1000, .combine=rbind) %dopar% {
      rdirichlet(1,c(3,4,2,0.5,0.5)) %>% PiToTheta(.) %>% t(.)
}

# marginal distribution
layout(matrix(1:4,nrow=2))
for (i in 1:4)
{hist(theta.prior[,i],breaks=25)}
layout(matrix(1,nrow=1))

# covariance
cor(theta.prior)
pairs(theta.prior)

##---------------------------------------------------------------------##
library(MVB)
n <- 1000
p <- 5
kk <- 2
tt <- NULL
alter <- 1
for (i in 1:kk) {
      vec <- rep(0, p)
      vec[i] <- alter
      alter <- alter * (-1)
      tt <- cbind(tt, vec)
}
tt <- 1.5 * tt
tt <- cbind(tt, c(rep(0, p - 1), 1))
x <- matrix(rnorm(n * p, 0, 4), n, p)
res <- mvb.simu(tt, x, K = kk, rep(.5, 2))
fitMVB <- mvbfit(x, res$response, output = 1)



