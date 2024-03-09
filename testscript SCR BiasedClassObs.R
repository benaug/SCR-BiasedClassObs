library(nimble)
library(coda)
source("sim.SCR.BCO.R")
source("init.SCR.BCO.R")
source("NimbleModel BCO.R")
source("NimbleFunctions BCO.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
n.class <- 2 #how many classes
pi.class <- c(0.45,0.55) #class frequencies
p.obs.class <- c(0.5,0.9) #class observation probabilities|capture
p0 <- c(0.05,0.3)
sigma <- c(1,0.5)
if(!(length(pi.class)==n.class&length(p0)==n.class&length(sigma)==n.class))stop("pi.class, p0, and sigma must be of length n.class")

K <- 5
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- as.matrix(expand.grid(3:11,3:11)) #trapping array

D <- 0.6 #see lambda below for what this D implies given state space area
area <- get.area(X,buff)
lambda <- D*area #expected N
lambda

data <- sim.SCR.BCO(D=D,p0=p0,sigma=sigma,pi.class=pi.class,p.obs.class=p.obs.class,
                    K=K,X=X,buff=buff)
#total individual caps
rowSums(data$y)
#total individual class observations
rowSums(data$y.obs.class)
#total individual caps by class
rowSums(data$y[data$class.obs==1,,])
rowSums(data$y[data$class.obs==2,,])

#Data augmentation level for each session
M <- 175
#initialize data structures
nimbuild <- init.SCR.BCO(data,M)

z.data <- rep(NA,M)
z.data[1:nimbuild$n] <- 1

J <- nimbuild$J
K1D <- rep(K,J) #trap operation

#get class data
fixed.class <- which(rowSums(data$y.obs.class)>0)
class.data <- rep(NA,M)
class.data[fixed.class] <- data$class.obs[fixed.class]

#inits for nimble
N.init <- sum(nimbuild$z.init)
class.init <- sample(1:n.class,M,replace=TRUE)
class.init[fixed.class] <- class.data[fixed.class]

Niminits <- list(z=nimbuild$z.init,s=nimbuild$s.init,classes=class.init,
                 p0=c(0.1,0.1),sigma=c(1,1),N=N.init,D=0.5,p.obs.class=c(0.5,0.5),pi.class=c(0.5,0.5),alpha=c(1,1))

#constants for Nimble
constants <- list(M=M,J=J,K1D=K1D,xlim=nimbuild$xlim,ylim=nimbuild$ylim,area=area,
                  n.class=n.class,n=nimbuild$n)

#supply data to nimble
y.obs.class2D <- apply(data$y.obs.class,c(1,2),sum)
#note, using "classes" instead of "class" in nimble. Using the latter wouldn't compile for some reason.
Nimdata <- list(y=nimbuild$y,y.obs.class=y.obs.class2D,z=z.data,X=nimbuild$X,classes=class.data)

# set parameters to monitor
parameters<-c('lambda','p0','sigma','N','D',"pi.class","p.obs.class")

nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE) 

###*required* sampler replacement
z.ups <- round(M*0.25) # how many z proposals per iteration per session?
conf$removeSampler("N")

#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames(paste("y[1:",M,",1:",J,"]"))
pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames("N")
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,y.nodes,pd.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(inds.detected=1:nimbuild$n,z.ups=z.ups,J=J,M=M,
                                                 xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),
                silent = TRUE)


#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1.
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim),silent = TRUE)
}

#adding a block sampler for p0 and sigma by class. Keeping independent samplers
# conf$removeSampler(c("p0","sigma"))
for(c in 1:n.class){
  conf$addSampler(target = c(paste("p0[",c,"]", sep=""),paste("sigma[",c,"]", sep="")),
                  type = 'RW_block',control=list(),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$N #realized N target
