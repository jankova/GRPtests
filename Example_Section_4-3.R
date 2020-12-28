gen_real_x <- function(n, p) {
  X <- rnorm(n*p); dim(X) <- c(n, p)
  rel_vars <- var_order[1:p]
  X <- X %*% sqrt_cov_x_trans[rel_vars, rel_vars]
  X <- pnorm(X)
  for (j in 1:p) {
    X[, j] <- quantile(x_orig[, rel_vars[j]], X[, j])
  }
  return(X)
}

library(glmnet)
library(GRPtests)

## LUNG CANCER DATASET gds2771 from NCBI database.
## n = 192, p = 22000; we reduce the number of predictors to p = 200

## LOAD THE DATASET

source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(Biobase)
library(GEOquery)
gds2771 <- getGEO('GDS2771', destdir=".")
dat <- gds2771@dataTable@table

## DEFINE Y TO BE ZERO-ONE VECTOR

y <- rep(0, nrow(gds2771@dataTable@columns))
for(i in 1:nrow(gds2771@dataTable@columns)){
  if(grepl("NOT", gds2771@dataTable@columns[i,3], fixed = TRUE) ||
     grepl("suspect", gds2771@dataTable@columns[i,3], fixed = TRUE)){
    y[i] <- 1 # 
  }
}

## take only rows without missing observations
dat <- dat[complete.cases(dat), ]
dat <- t(dat[,-c(1,2)])

###### GENERATE ARTIFICIAL Y 

## select 500 variables with the highest variances
vars <- apply(dat, 2, var)
max.vars500 <- dat[, sort(vars, decreasing = TRUE, index.return = TRUE)$ix <= 500]
x_orig <- scale(max.vars500)

n <- nrow(x_orig)
p <- ncol(x_orig)

# Transformation to normality
cov_x_trans <- cov(qnorm(apply(x_orig, 2, rank) / (n+1)))
sqrt_cov_x_trans <- chol(cov_x_trans, pivot=TRUE)
pivot <- attr(sqrt_cov_x_trans, "pivot")
oo <- order(pivot)
sqrt_cov_x_trans <- sqrt_cov_x_trans[, oo]
var_order <- order(apply(x_orig, 2, var), decreasing=TRUE)

## generate new Gaussian data
Xn <- gen_real_x(800, p)

## Fit Lasso to estimate the support set

glmnet.fit <- cv.glmnet(x_orig, y, family = "binomial", alpha = 0.5)  
summary(glmnet.fit)
coefs <- as.vector(coef(glmnet.fit))
S.hat <- which(coefs != 0)

## We consider 3 settings as in the paper: 
## g(u) = 0, g(u) = u_{j_1}^2 + u_{j_2}^2, g(u) = u_{j_1}*u_{j_2} + u_{j_3}*u_{j_4}

## NOTE: THE SETTING MUST CHOSEN MANUALLY FROM THE FOLLOWING THREE OPTIONS:

## Setting 1: g(u) = 0 
g <- 0

## Setting 2: g(u) = u_{j_1}^2 + u_{j_2}^2
set.seed(1)
ids <- sample(S.hat, 2)
newvar <- Xn[,ids[1]]^2 + Xn[,ids[2]]^2 
g <- sqrt(3) * scale(newvar)

## Setting 3: g(u) = u_{j_1}*u_{j_2} + u_{j_3}*u_{j_4}
set.seed(1)
ids2 <- sample(S.hat, 4)
newvar <- Xn[, ids2[1]] * Xn[, ids2[2]] + Xn[, ids2[3]] * Xn[, ids2[4]]
g <- sqrt(3) * scale(newvar)

#######################################

# calculate probabilities
pr = 1 / (1 + exp( - (Xn %*% coefs[-1] + coefs[1] + g )))

iter <- 100
pv <- rep(0, iter)

set.seed(1)
for(i in 1:iter){
  
  print(i)
  yn <- rbinom(nrow(Xn), 1, pr)

  test.obj <- GRPtest(Xn, yn, fam = "binomial", nsplits = 1) 
  pv[i] <- test.obj
  print(test.obj)

}

mean(pv < 0.05)

## END OF CODE. ##

