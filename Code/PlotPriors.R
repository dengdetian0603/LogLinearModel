K = 5
Sigma = 1.6
Alpha = 0.8
x = rnorm(1000*K, 0, Sigma)
mu = matrix(exp(x)/(1+exp(x)),ncol=K)
pi0 = apply(mu, 1, function(mu) runif(1,0,1-max(mu)))
pi1toK = sapply(pi0, function(x) {tmp=rStickBreak(K-1, Alpha); 
                                    tmp2=c(tmp,1-sum(tmp)); 
                                    return(tmp2*(1-x))})
prior = cbind(mu,pi0,t(pi1toK))
apply(prior[,-(1:K)],2,mean)


layout(matrix(1:2,ncol=1))
par(mar=c(2,4,1,1))
plot(density(prior[,1],adjust=1.8), xlab="value of parameter", xlim=c(0,1),main="Posterior Density", ylim=c(0,7.5))
lines(density(prior[,2],adjust=1.8),col=2)
lines(density(prior[,3],adjust=1.8),col=4)
legend("topright",legend=c(paste0("Mu1 : ",round(mean(prior[,1]),3)),
                           paste0("Mu2 : ",round(mean(prior[,2]),3)),
                           paste0("Mu3 : ",round(mean(prior[,3]),3))),
       lty=1, col=c(1,2,4))
plot(density(prior[,K+1],adjust=1.8), xlab="value of parameter", xlim=c(0,1), ylim=c(0,6), main="Posterior Density")
lines(density(prior[,K+2],adjust=1.8),col=2)
lines(density(prior[,K+3],adjust=1.8),col=4)
lines(density(prior[,K+4],adjust=1.8),col=3)
legend("topright",legend=c(paste0("Pi0 : ",round(mean(prior[,K+1]),3)),
                           paste0("Pi1 : ",round(mean(prior[,K+2]),3)),
                           paste0("Pi2 : ",round(mean(prior[,K+3]),3)),
                           paste0("Pi3 : ",round(mean(prior[,K+4]),3))),
       lty=1, col=c(1,2,4,3))
layout(matrix(1))

round(apply(prior,2,mean),3)

