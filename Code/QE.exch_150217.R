#  Quadratic Exponential Family (QE)
## Exchangeable parameterization
obj = function(K,pi0,pi1,gamma)
{
      value = 0
      for (s in 2:K)
      {
            value = value + pi0*choose(K,s)*exp(s*log(pi1/pi0/K)+choose(s,2)*gamma)
      }
      value + pi0 + pi1 -1
}

obj(K=4,pi0=0.4,pi1=0.3,gamma=10)


probs = function(K,pi0,pi1,gamma)
{
      Pi = c(pi0,pi1)
      for (s in 2:K)
      {
            Pi[s+1] = pi0*choose(K,s)*exp(s*log(pi1/pi0/K)+choose(s,2)*gamma)
      }
      Pi
}

ga = uniroot(f=obj,interval=c(-200,200),K=10,pi0=0.3,pi1=0.4)
Pi = probs(K=10,pi0=0.3,pi1=0.4,gamma=ga$root)
barplot(Pi)

plot(0:10,log(Pi))