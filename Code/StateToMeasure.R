library(doMC)
library(foreach)

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

