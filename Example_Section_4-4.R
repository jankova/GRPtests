## Group testing in logistic regression: Comparison of GRPgrouptest, globaltest and de-biased lasso

gen.data <- function(n, p, fam, sig){
  s <- 4
  beta <- 1*c(rep(1,min(p,s)),rep(0,p-min(p,s)))    # parameters
  beta[5] <- sig
  x <- rho^c(0:(p-1))
  sigma0 <- toeplitz(x)
  
  X <- mvrnorm(n, rep(0,p), sigma0)*1
  
  z = crossprod(t(X), beta) 
  
  if(fam == "gaussian"){
    y = z + rnorm(n)
  }else{
    pr = 1/(1 + exp(-z)) 
    y = rbinom(n, 1, pr)      
  }
  
  return(list(X = X, y = y))
}

plots <- function(){
  plot(ecdf(pvals.list[[1]]), col = "red", xlim = c(0,1), ylim=c(0,1))
  lines(ecdf(pvals.list[[2]]), col="blue")  
  #lines(ecdf(pvals.list[[3]]), col="black")
  legend("bottomright", c("GRP-test", "Globaltest", "De-biased Lasso"), 
         col = c("red", "blue", "black"), pch = c(20,20,20), lty = c(1,1,1), lwd = c(2,2,2), bty = 'n')
  abline(0,1)
  
}

all.tests <- function(sig){
  iter <- 10
  pvals.grouptest <- pvals.globaltest <- pvals.lasso <- rep(0, iter)
  for(i in 1:iter){
    
    print(paste(" NUMITER: ",i,".. "))
    gd <- gen.data(n, p, fam, sig)
    X <- gd$X
    y <- gd$y
    
    # globaltest
    gtest <- gt(y ~ X[, c(1:p)[-G]], y ~ X[,c(1:p)], model = "logistic")
    pvals.globaltest[i] <- gtest@result[1]
    
    # GRPgrouptest
    pvals.grouptest[i] <- GRPgrouptest(X, y, fam = "binomial", G = G, B = 100, penalize = TRUE)
    
    # de-biased lasso grouptest (NOTE THAT THIS TAKES A LOT OF COMPUTATIONAL TIME!)
    #fit.lasso <- lasso.proj(X, y, family = fam, suppress.grouptesting = FALSE)
    #pvals.lasso[i] <- fit.lasso$groupTest(G)
    
  }
  
  list(pvals.grouptest, pvals.globaltest, pvals.lasso)
}


library(GRPtests)
library(globaltest)

## test group G
G <- c(5:100)
fam <- "binomial"
rho <- 0.6

set.seed(1)

## Compare GRP-test, global and lasso, n = 500, p = 100

n = 500
p = 100
## UNDER NULL
sig = 0
pvals.list <- all.tests(sig)
plots()

## UNDER ALTERNATIVES 
sig = 1
pvals.list <- all.tests(sig)
plots()

## Compare GRP-test and globaltest, n = 500, p = 800


n = 600
p = 800
G = 5:p
## UNDER NULL
sig <- 0
pvals.list <- all.tests(sig)
plots()

## UNDER ALTERNATIVES
sigs <- seq(0, 1.4, by = 0.2)
powers <- matrix(0, nrow = length(sigs), ncol = 2)

for(j in 1:length(sigs)){
  sig <- sigs[j]
  pvals.list <- all.tests(sig)
  powers[j, 1] <- mean( pvals.list[[1]] < 0.05) 
  powers[j, 2] <- mean( pvals.list[[2]] < 0.05) 
}

plot(sigs, powers[,1], col ="red", ylim = c(0,1), xlim = c(0, 1.4), type="l", lwd = 2,
     ylab= "Probability of rejection", xlab = "")
lines(sigs, powers[,2], col ="blue", lwd = 2)

