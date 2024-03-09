NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  D ~ dunif(0,100) #expected D
  for(c in 1:n.class){
    p0[c] ~ dunif(0,1)
    sigma[c] ~ dunif(0,100)
    p.obs.class[c] ~ dunif(0,1)
    log(alpha[c]) ~ dnorm(0,sd=2)
  }
  pi.class[1:n.class] ~ ddirch(alpha[1:n.class])
  #--------------------------------------------------------------
  lambda <- D*area #expected N
  N ~ dpois(lambda) #realized N
  for(i in 1:M){
    classes[i] ~ dcat(pi.class[1:n.class])
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma[classes[i]], p0=p0[classes[i]], z=z[i])
    y[i,1:J] ~ dBernoulliVector(pd=pd[i,1:J],K1D=K1D[1:J],z=z[i]) #vectorized obs mod
  }
  #Class observation process
  for(i in 1:n){
    for(j in 1:J){
      y.obs.class[i,j] ~ dbinom(size=y[i,j],prob=p.obs.class[classes[i]])
    }
  }
})