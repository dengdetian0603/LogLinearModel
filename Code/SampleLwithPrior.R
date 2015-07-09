library(limSolve)
library(doMC)
library(foreach)
source("~/Documents/Johns Hopkins SPH/Research/S.Zeger/LogLinearModel/Code/ToolBox_PERCH.R")

registerDoMC(4) 

K = 5 # 
mu.prior = c(0.79,0.79,0.79,0.79,0.79)
pi.prior = c(0.2,0.5,0.298998,0.001,0.000001,0.000001)


#K = 2 # when K < 3, only single solution
#mu.prior = c(0.4,0.3)
#pi.prior = c(0.4,0.5,0.1)

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

## Given Pi only
Phis = xsample( E = PiMat, F = pi.prior[-1]/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                iter = 6000, output = 100, type = "cda")
mu.sample = MuMat%*%t(Phis$X[1:100,])*pi.prior[1]
boxplot(t(mu.sample))

## Given Mu and Pi0 only
# need additional constrain that pis sum up to 1
Phis = xsample( E = rbind(rep(1,2^K-1),MuMat), F = c(1-pi.prior[1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
                iter = 6000, output = 100, type = "mirror")
pi.sample = PiMat%*%t(Phis$X[1:100,])*pi.prior[1]
#round(pi.sample,3)
boxplot(t(pi.sample))

## Given both Mu and Pi
mu.prior = mu.sample[,1]
Phis = xsample( E = rbind(PiMat,MuMat), F = c(pi.prior[-1],mu.prior)/pi.prior[1], G = diag(rep(1,J1)), H = rep(0,J1),
               iter = 300, output = 100, type = "cda")
Phis$X

L.all

i=10
cell.prob = c(1, Phis$X[i,])*pi.prior[1]
rmultinom(1,100,cell.prob)

apply(L.all,1,function(x) paste(as.character(x),collapse=""))

rDataWithPiMu.LL(5, mu.prior, pi.prior,1)

## solve for qe prameters
library(nloptr)

# Objective function
eval_f = function(x, xx, B, Lmat, b) 
{
      list("objective" = sum((x-xx)^2),
           "gradient" = 2*(x-xx))
}

# non-linear equality constraints
eval_g_eq = function(x, B, Lmat, b,xx)
{
      constr = B%*%exp(Lmat%*%x) - b
      grad = B%*%diag(exp(Lmat%*%x)[,1])%*%Lmat
      return( list( "constraints"=constr, "jacobian"=grad ) )
}

# initial values
npar = (K*(K+1)/2)
x0 = Phis$X[2,1:npar]
# lower and upper bounds of control
lb = rep(-Inf,npar)
ub = rep(Inf,npar)
local_opts = list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-5 )
opts = list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-5,
              "maxeval" = 20000,
              "local_opts" = local_opts )
res = nloptr( x0=x0, eval_f=eval_f,lb=lb, ub=ub,
              eval_g_eq=eval_g_eq, opts=opts,
              xx=x0, B=rbind(PiMat,MuMat), b=c(pi.prior[-1],mu.prior)/pi.prior[1], Lmat = L.all.exp[,1:npar])
res


