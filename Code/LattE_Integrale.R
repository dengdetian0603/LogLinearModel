source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/ToolBox_PERCH.R")
source("/home/bst/student/ddeng/ThesisTopic/LogLinearModel/GitRepo/LogLinearModel/Code/Sparse_Models.R")
#library(MVB)
#library(limSolve)
library(doMC)
library(foreach)
#library(coda)
library(MASS)
#library(editrules)

num_core = min(detectCores(), 12)
registerDoMC(num_core)

load("SimDat_K5n300.RData")

Lmat = GenMatrices.Sparse(K=5, Smax=3)
MuMat = Lmat$MuMat
PiMat = Lmat$PiMat
J1 = Lmat$J1


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

tryCatch(xsample( E = rbind(pars$PiMat,pars$MuMat), F = c(pis[-1],mu)/pis[1], G = diag(rep(1,pars$J1)), H = rep(0,pars$J1),
                      iter = n, burnin = burnin, type = method, test=FALSE), error = function(e) "Incompatible constraints.")


# ------------------------------------------------------------------------------------------ #
# call LattE function integrate from R
mu.tmp = c(0.258,0.157,0.165,0.288,0.541)
pi.tmp = c(0.103,0.513,0.256,0.128)


