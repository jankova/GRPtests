## High-dimensional logistic regression: Goodness-of-fit testing
library(MASS)
library(glmnet)
library(randomForest)
library(globaltest)

gen.data <- function(n, p, fam, quad){
  s <- 5
  beta <- 1*c(rep(1,min(p,s)),rep(0,p-min(p,s)))    # parameters
  x <- rho^c(0:(p-1))
  sigma0 <- toeplitz(x)
  
  X <- mvrnorm(n, rep(0,p), sigma0)*1#matrix(runif(n*p),n,p)#
  
  if(quad == TRUE){
    f =  sig * (X[,1]^2) #*X[,3])#^2 -1)#-sigma0[1,2])#*(X[,1]^2-1) #*X[,5]#^2
  }else{
    f = sig * X[,1] * X[,2]
  }
  
  z = f + crossprod(t(X), beta) 
  
  if(fam == "gaussian"){
    y = z +rnorm(n)#rlaplace(n, m=0, s=5)#1*rnorm(n)
  }else{
    pr = 1/(1 + exp(-z)) 
    # pr <- pnorm(z)                    # misspecified link: probit
    y = rbinom(n, 1, pr)      
  }
  
  return(list(X = X, y = y))
}


n = 900
p = 900
it10 <- 200
fam <- "binomial"

# QUADRATIC
sigs <- seq(0, 1.5, by = 0.5)#
quad <- TRUE

# Interaction
sigs <- seq(0, 3, by = 1)
quad <- FALSE

pv <- array(0,dim = c(3,length(sigs),it10))
pvbench <- array(0,dim = c(3,length(sigs),it10))

rhos<-c(0.4,0.6,0.8)
irhos <- c(1,2,3)


set.seed(1)
for(j in 1:length(rhos)){
  rho <- rhos[j]
  #n <- ns[j]
  ir <- irhos[j]
  print("")
  print(paste("rho:",rhos[j]))
  for(k in 1:length(sigs)){
    sig <- sigs[k]
    print(paste("sig:",sigs[k]))
    set.seed(2)
    for(i in 1:it10){
      cat(paste(i,".."))
      gd <- gen.data(n, p, fam, quad)
      X <- gd$X
      y <- gd$y
      
      obj.pval <- GRPtest(X, y, fam = "binomial", nsplits = 1)
      obj.pval.bench <- GRPtest(X[, 1:5], y, fam = "binomial", nsplits = 1)
      
      pv[ir, k, i] <- obj.pval #test.obj$pval
      pvbench[ir, k, i] <- obj.pval.bench#test.obj.bench$pval
      cat(paste(round(obj.pval,2),"=pv"))
      cat(paste(round(obj.pval.bench,2),"=bench"))
    }
  }
}

pvbench_rho1 <- apply(pvbench[1,,], 1, function(x) mean(x < 0.05))
pvbench_rho2 <- apply(pvbench[2,,], 1, function(x) mean(x < 0.05))
pvbench_rho3 <- apply(pvbench[3,,], 1, function(x) mean(x < 0.05))

p_rho1 <- apply(pv[1,,], 1, function(x) mean(x < 0.05))
p_rho2 <- apply(pv[2,,], 1, function(x) mean(x < 0.05))
p_rho3 <- apply(pv[3,,], 1, function(x) mean(x < 0.05))

library(xtable)
titl <- c("sigma = 0", "sigma = 0.5", "sigma = 1", "sigma = 1.5")
xtable(rbind(titl, pvbench_rho1, pvbench_rho2, pvbench_rho3),
       align = NULL, digits = NULL,
       display = NULL, auto = FALSE)

xtable(rbind(titl, p_rho1, p_rho2, p_rho3),
       align = NULL, digits = NULL,
       display = NULL, auto = FALSE)


par(mfrow=c(length(rhos),length(sigs)))
for(j in 1:length(rhos)){
  for(k in 1:length(sigs)){
    plot(ecdf(pv[j, k, ]),xlim=c(0,1),ylim=c(0,1))
    abline(0,1)
  }
}

### END OF CODE. ###