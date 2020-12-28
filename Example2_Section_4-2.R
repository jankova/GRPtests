## High-dimensional logistic regression: Goodness-of-fit testing
library(MASS)
library(glmnet)
library(randomForest)
library(globaltest)
library(GRPtests)
library(mvnfast)

## AUXILIARY FUNCTIONS

gen.data <- function(n, p, fam, quad, sigma0, beta){

  #X <- mvrnorm(n, rep(0,p), sigma0)*1
  X <- rmvn(n, rep(0,p), sigma0)*1
  
  if(quad == TRUE){
    f =  sig * X[,1]^2
  }else{
    f = sig * X[,1] * X[,2]
  }
  
  z = f + crossprod(t(X), beta) 
  
  if(fam == "gaussian"){
    y = z +rnorm(n)
  }else{
    pr = 1/(1 + exp(-z)) 
    y = rbinom(n, 1, pr)      
  }
  
  return(list(X = X, y = y))
}



nonlin_test <- function(X, y, fam, nsplits, RP_function, penalize){
  ## Computes p-values for goodness-of-fit test from "nsplits" splits of data
  ## and aggregates them.
  
  pvals <- rep(0, nsplits)
  
  for(i in 1:nsplits){
    # Compute p-value from a single split
    pvals[i] <- single_split_pval(X, y, fam, RP_function, penalize)
  }
  print(pvals)
  if(nsplits > 1){
    pval <- min(2*median(pvals), 1)
  }else{
    pval <- pvals[1]
  }
  return(list(pval = pval, pvals = pvals))
  
}

single_split_pval <- function(X, y, fam, RP_function, penalize){
  ## Computes a p-value for goodness-of-fit test from a single split of data.
  # Split sample
  p <- ncol(X)
  n <- nrow(X)
  ind <- sample.int(n = nrow(X), size = floor(.5*nrow(X)), replace = FALSE)
  
  XA <- as.matrix(X[ind, ])
  XB <- as.matrix(X[-ind, ])
  yA <- y[ind]
  yB <- y[-ind]
  
  # Fit GLMNET on Part A
  
  partA <- glmfit(XA, yA, fam, penalize)
  partB <- glmfit(XB, yB, fam, penalize)
  
  # Compute residuals resA, resB
  resA <- partA$res/partA$D
  resB <- partB$res/partB$D
  
  beta.hat <- partA$beta #### CHANGE
  hatS <- which(beta.hat[-1] != 0)
  S100 <- partA$S100
  
  # Compute D
  mix <- XB%*%beta.hat[-1] + beta.hat[1]  # predict fitA on XB
  pars <- GLM_parameters(mix, fam)
  D <- as.vector(sqrt(pars$var_y))
  
  # Apply RP-function to residuals
  
  if (fam == "binomial" && sum(is.nan(resA) != 0)){
    stop("Residuals from glmnet contain NaN. One binomial class might have too few observations.")
  }
  
  pred.rf <- RP_function(XA, partA$res, XB, S100)
  
  # Orthogonalize
  w <- orthogonalize(D, XB, pred.rf, hatS)
  
  w <- w / as.numeric(sqrt(crossprod(w)))
  
  # Test statistic and p-value
  TS <- as.numeric(crossprod(w, resB))
  pval <- 1 - pnorm(TS)
  
  pval
}

orthogonalize <- function(D, XB, pred.rf, hatS){
  ## Uses sqrt_lasso to approximately orthogonalize direction pred.rf
  ## with respect to variables in XB.
  ## Exact orthognalization is done on variables contained in hatS.
  
  ## Fit sqrt-LASSO to estimate w
  Dw <- D
  p <- ncol(XB)
  n <- nrow(XB)
  
  # Define indices 'indi' where we orthognalize exactly
  indi <- 1*rep(1,p)
  indi[hatS] <- 0
  if(sum(indi) == 0) indi[1] <- 1 # if all penalty factors are zero, put penalty on the first entry to avoid error
  
  # Orthogonalize
  fit <- glmnet(Dw*XB, Dw*pred.rf, penalty.factor = indi, nlambda = 40)
  W <- Dw*pred.rf - Dw*XB%*%fit$beta
  lambda.sqrt.lasso <- sqrt(2*log(p)/n)
  index <- which.min(abs(fit$lambda - lambda.sqrt.lasso*apply(W, 2, function(x) sqrt(sum(x^2))/sqrt(n) )))
  beta.sqrt <- fit$beta[,index]
  
  w <- Dw*pred.rf - Dw*XB%*%beta.sqrt
  w
  
  #max(abs(t(XB)%*%diag(Dw^2)%*%w))/sqrt(crossprod(w))
  #sqrt(log(p))
}


glmfit <- function(XA, yA, fam, penalize){
  ## Fits glmnet to data (XA, yA) and output several values from the fit
  
  if(fam == "binomial"){
    y.trans <- as.factor(yA)
  }else{
    y.trans <- yA
  }
  
  if(penalize == FALSE){
    fit1 <- glmnet(XA, y.trans, family = fam, lambda = 0)
    beta.hat <- coef(fit1)
    mix <- XA %*% beta.hat[-1] + beta.hat[1]
    S100 <- 1:ncol(XA)
  }else{
    fit1 <- cv.glmnet(x = XA, y = y.trans, family = fam)
    beta.hat <- coef(fit1, s = "lambda.min")
    mix <- XA %*% beta.hat[-1] + beta.hat[1]
    
    temp <- c(1:length(fit1$glmnet.fit$df))[fit1$glmnet.fit$df < 100]
    
    index.lam <- temp[length(temp)]
    
    lam100 <- fit1$glmnet.fit$lambda[index.lam]
    beta.hat100 <- coef(fit1, s = lam100)
    S100 <- which(beta.hat100[-1] != 0)
  }
  
  
  pars <- GLM_parameters(mix, fam)
  mean_y <- pars$mean_y
  var_y <- pars$var_y
  D <- sqrt(var_y)
  
  resA <- yA - mean_y
  
  list(D = as.vector(D), res = as.vector(resA), beta = as.vector(beta.hat), 
       mu = as.vector(mean_y), S100 = S100)
  
}

RP_randomForest <- function(XA, resA, XB, S100){
  ## Fits a random forest of resA on XA and predicts it on XB
  
  #rf <- randomForest(XA[, vars], as.vector(resA))
  #pred.rf <- predict(rf, newdata = XB[, vars])
  
  vars <- 1:ncol(XA)
  num.vars <- length(vars)
  data.res <- data.frame(rbind(XA[, vars], XB[, vars]), c(resA, resA))
  rf <- ranger(y = data.res[1:nrow(XA),num.vars + 1], 
               x = data.res[1:nrow(XA),1:num.vars], data = data.res[1:nrow(XA),])
  pred.rf <- predict(rf, data = data.res[(nrow(XA)+1):(nrow(XA)+nrow(XB)),])$predictions
  
  pred.rf
}


GLM_parameters <- function(mix, fam){
  ## Computes predicted mean and variance of a GLM
  if(fam == "gaussian"){
    mean_y <- mix
    var_y <- rep(1, nrow(mix))
  }else if(fam == "binomial"){
    mean_y <- exp(mix) / (1 + exp(mix))
    var_y <- mean_y * (1 - mean_y)
  }else if(fam == "poisson"){
    mean_y = exp(mix)
    var_y = exp(mix)
  }
  list(mean_y = as.vector(mean_y), var_y = as.vector(var_y))
}

GRPtest <- function(X, y, fam = c("gaussian", "binomial", "poisson"),
                    RP_function = NULL,
                    nsplits = 5L, penalize = ifelse(p >= floor(n/1000), TRUE, FALSE), output_all = FALSE){
  
  if (!is.matrix(X)){
    stop("X should be a matrix with at least one column.")
  }
  np <- dim(X)
  if (is.null(np) | (np[2] < 1L)){
    stop("X should be a matrix with at least one column.")
  }
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  y <- as.numeric(y)
  if (length(y) != n){
    stop("y must have nrow(X) components.")
  }
  
  if ((p >= n-1) && (penalize == FALSE)){
    stop("When penalize=FALSE we must have ncol(X) < nrow(X)-1; try setting penalize=TRUE.")
  }
  
  if(fam != "gaussian" && fam != "binomial" && fam != "poisson"){
    stop("family must be gaussian, binomial or poisson; try changing parameter fam.")
  }
  
  if (is.null(RP_function)) {
    RP_function <- RP_randomForest
  }
  output <- nonlin_test(X, y, fam, nsplits, RP_function, penalize)
  pval_computed <- output$pval
  if (output_all == FALSE){
    return(pval_computed)
  }else{
    return(output)
  }
}


############################
### MAIN CODE ###
############################

n = 2000
p = 3000
it10 <- 100
fam <- "binomial"

# QUADRATIC
sig_quad <- 1.5

# Interaction
sig_int <- 3

s <- 5
beta <- 1*c(rep(1,min(p,s)),rep(0,p-min(p,s)))    # parameters
x <- rho^c(0:(p-1))
sigma0 <- toeplitz(x)

rho <- 0.6

# no misspecification
pv <- pvbench <- pvbenchP <- rep(0, it10)

# nonlinearity: quadratic term 
pv_quad <- pvbench_quad <- rep(0, it10)

# nonlinearity: interaction
pv_int <- pvbench_int <- rep(0, it10)


start <- Sys.time()
set.seed(1)
for(i in 1:it10){
  print(paste(i,".."))

  sig <- 0
  gd <- gen.data(n, p, fam, quad = TRUE, sigma0, beta)
  X <- gd$X
  y <- gd$y
  obj.pval <- GRPtest(X, y, fam = "binomial", nsplits = 1, penalize = TRUE)
  obj.pval.bench <- GRPtest(X[, 1:5], y, fam = "binomial", nsplits = 1, penalize = FALSE)
  obj.pval.benchP <- GRPtest(X[, 1:5], y, fam = "binomial", nsplits = 1, penalize = TRUE)
  
  pv[i] <- round(obj.pval,20)
  pvbench[i] <- round(obj.pval.bench,20)
  pvbenchP[i] <- round(obj.pval.benchP,20)
  
  sig <- 1.5
  gd <- gen.data(n, p, fam, quad = TRUE, sigma0, beta)
  X <- gd$X
  y <- gd$y
  obj.pval <- GRPtest(X, y, fam = "binomial", nsplits = 1)
  obj.pval.bench <- GRPtest(X[, 1:5], y, fam = "binomial", nsplits = 1)
  
  pv_quad[i] <- round(obj.pval,20)
  pvbench_quad[i] <- round(obj.pval.bench,20)
  
  
  sig <- 2
  gd <- gen.data(n, p, fam, quad = FALSE, sigma0, beta)
  X <- gd$X
  y <- gd$y
  obj.pval <- GRPtest(X, y, fam = "binomial", nsplits = 1)
  obj.pval.bench <- GRPtest(X[, 1:5], y, fam = "binomial", nsplits = 1)
  
  pv_int[i] <- round(obj.pval,20)
  pvbench_int[i] <- round(obj.pval.bench,20)
  
}
Sys.time()-start



plot(ecdf(pv[1:i]),xlim=c(0,1), col="red", xlab="p-values", ylab = "Ecdf of p-values",
     main = "")
lines(ecdf(pvbench[1:i]),xlim=c(0,1))
legend("bottomright", c("GRP-test", "Benchmark"), col = c("red", "black"),
       pch=c(20,20), lwd = 1, bty ='n')
abline(0,1)

plot(ecdf(pv_quad[1:i]),xlim=c(0,1), col="red", xlab="p-values", ylab = "Ecdf of p-values",
     main = "")
lines(ecdf(pvbench_quad[1:i]),xlim=c(0,1))
legend("bottomright", c("GRP-test", "Benchmark"), col = c("red", "black"),
       pch=c(20,20), lwd = 1, bty ='n')
abline(0,1)

plot(ecdf(pv_int),xlim=c(0,1), col="red", xlab="p-values", ylab = "Ecdf of p-values",
     main = "")
lines(ecdf(pvbench_int),xlim=c(0,1))
legend("bottomright", c("GRP-test", "Benchmark"), col = c("red", "black"),
       pch=c(20,20), lwd = 1, bty ='n')
abline(0,1)


### END OF CODE. ###

#pv_all <- list(pv, pvbench, pv_quad, pvbench_quad, pv_int, pvbench_int)

#saveRDS(list(n, p, rho, sigma0, sig, beta,pv_all), file = "n2000p3000_sig0rho06new2.rds")
