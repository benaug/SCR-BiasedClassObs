init.SCR.BCO <-   function(data=data,M=M){
  J <- nrow(data$X)
  y <- matrix(0,M,J)
  n <- nrow(data$y)
  y[1:n,] <- apply(data$y,c(1,2),sum)
  s.init <- cbind(runif(M,data$xlim[1],data$xlim[2]), runif(M,data$ylim[1],data$ylim[2])) #assign random locations
  idx <- which(rowSums(y)>0) #switch for those actually caught
  for(i in idx){
    trps <- matrix(X[y[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  z.init <- 1*(rowSums(y)>0)
  return(list(y=y,X=data$X,xlim=data$xlim,ylim=data$ylim,K=data$K,J=J,n=n,s.init=s.init,z.init=z.init))
}
