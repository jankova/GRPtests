
nonlin_test <- function(X, y, fam, nsplits, RP_function, penalize){
## Computes p-values for goodness-of-fit test from "nsplits" splits of data
## and aggregates them.

  pvals <- rep(0, nsplits)

  for(i in 1:nsplits){
    # Compute p-value from a single split
    pvals[i] <- single_split_pval(X, y, fam, RP_function, penalize)
  }

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
  resB <- partB$res/partA$D  

  beta.hat <- partA$beta
  hatS <- which(beta.hat[-1] != 0)

  # Compute D
  mix <- XB%*%beta.hat[-1] + beta.hat[1]  # predict fitA on XB
  pars <- GLM_parameters(mix, fam)
  D <- as.vector(sqrt(pars$var_y))

  # Apply RP-function to residuals

  if (fam == "binomial" && sum(is.nan(resA) != 0)){
    stop("Residuals from glmnet contain NaN. One binomial class might have too few observations.")
  }

  pred.rf <- RP_function(XA, resA, XB)

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
  indi <- rep(1,p)
  indi[hatS] <- 0
  if(sum(indi) == 0) indi[1] <- 1 # if all penalty factors are zero, put penalty on the first entry to avoid error

  # Orthogonalize
  fit <- glmnet(Dw*XB, Dw*pred.rf, penalty.factor = indi, nlambda = 40)
  W <- Dw*(pred.rf - XB%*%fit$beta)
  lambda.sqrt.lasso <- sqrt(2*log(p)/n)
  index <- which.min(abs(fit$lambda - lambda.sqrt.lasso*apply(W, 2, function(x) sqrt(sum(x^2))/sqrt(n) )))
  beta.sqrt <- fit$beta[,index]

  w <- Dw*(pred.rf - XB%*%beta.sqrt)
  w
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
  }else{
    fit1 <- cv.glmnet(x = XA, y = y.trans, family = fam)
    beta.hat <- coef(fit1, s = "lambda.min")
    mix <- XA %*% beta.hat[-1] + beta.hat[1]
  }

  pars <- GLM_parameters(mix, fam)
  mean_y <- pars$mean_y
  var_y <- pars$var_y
  D <- sqrt(var_y)

  resA <- yA - mean_y

  list(D = as.vector(D), res = as.vector(resA), beta = as.vector(beta.hat), mu = as.vector(mean_y))

}

RP_randomForest <- function(XA, resA, XB){
## Fits a random forest of resA on XA and predicts it on XB
  
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


### Function for group testing

group_test <- function(X, y, G, fam, B, penalize){
## Computes p-value for testing significance of predictors in G.

  n <- nrow(X)
  p <- ncol(X)
  W <- matrix(0, n, length(G))

  # Fit y on X[,-G]

  fit <- glmfit(X[,-G], y, fam, penalize)
  beta.hat <- fit$beta
  D <- fit$D
  mu <- fit$mu
  res <- ( y - mu ) / D

  # Compute test statistic
  test <- rep(0, length(G))

  for(j in 1:length(G)){
    k <- G[j]
    beta.sqrt <- RPtests::sqrt_lasso(as.matrix(D*X[,-G]), as.vector(D*X[,k]))
    w <- D*(X[,k] - X[,-G]%*%beta.sqrt)
    w <- w/as.numeric(sqrt(crossprod(w)))
    W[,j] <- w
    test[j] <- crossprod(w, res)
  }
  maxtest <- max(test)

  ## Multiplier bootstrap

  TB <- rep(0,B)
  for(i in 1:B){
    e <- rnorm(n)
    TB[i] <- apply((t(W*e)%*%res), 2, max)
  }

  pval <- (1 + sum(TB > maxtest))/(B + 1)
  return(pval)
}




