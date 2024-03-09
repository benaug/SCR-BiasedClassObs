e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
get.area = function (X, buff){
  area <- diff(range(X[,1])+c(-buff,buff))*diff(range(X[,2])+c(-buff,buff))
  return(area)
}
sim.SCR.BCO <-
  function(D=NA,p0=NA,sigma=NA,pi.class=NA,p.obs.class=NA,
           K=NA,X=X,buff=3){
    if(length(p0)!=length(sigma))stop("p0 and sigma must be of same length, one element per class")
    if(length(p.obs.class)!=length(sigma))stop("obsprob must be of length n.class")
    #simulate a population size
    area <- get.area(X,buff)
    lambda <- D*area #expected N
    N <- rpois(1,lambda) #realized N
    # simulate activity centers
    xlim <- range(X[,1])+c(-buff,buff)
    ylim <- range(X[,2])+c(-buff,buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    #simulate classes
    n.class <- length(p0)
    class <- sample(1:n.class,N,replace=TRUE,prob=pi.class)
    #get detection probs by class
    D <- e2dist(s,X)
    J <- nrow(X)
    pd <- p0[class]*exp(-D*D/(2*sigma[class]^2))
    
    # Capture individuals and observe their classes
    y <- array(0,dim=c(N,J,K))
    y.obs.class <- array(0,dim=c(N,J,K))
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          y[i,j,k] <- rbinom(1,1,pd[i,j])
          y.obs.class[i,j,k] <- rbinom(1,y[i,j,k],p.obs.class[class[i]])
        }
      }
    }

    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught <- which(apply(y,c(1),sum)>0)
    y.true <- y
    y <- y[caught,,]
    y.obs.class <- y.obs.class[caught,,]
    class.obs <- class[caught]
    if(K==1){
      y <- array(y,dim=c(dim(y),1))
      y.obs.class <- array(y.obs.class,dim=c(dim(y.obs.class),1))
    }
    n <- length(caught)
    out <- list(y=y,y.obs.class=y.obs.class,X=X,K=K,buff=buff,s=s,n=nrow(y),K=K,class=class,class.obs=class.obs,
              xlim=xlim,ylim=ylim,N=N)
    return(out)
  }
