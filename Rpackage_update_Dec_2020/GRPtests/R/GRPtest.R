#' Goodness-of-fit test for high-dimensional generalized linear models
#'
#' The function can test goodness-of-fit of a low- or high-dimensional
#' generalized linear model (GLM) by detecting the presence of nonlinearity in
#' the conditional mean function of y given X. Outputs a p-value.
#'
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param y Response vector.
#' @param fam Must be "gaussian", "binomial" or "poisson".
#' @param RP_function (optional) User specified function for residual prediction (see Details below).
#' @param penalize \code{TRUE} if penalization should be used when fitting the GLM models (see Details below).
#' @param nsplits Number of splits of the data set (see Details below).
#' @param output_all If \code{TRUE}, outputs all p-values from \code{nspilts} splits of the data.
#' @details
#'   This function tests if the conditional mean of \code{y} given \code{X} could be
#'   originating from a GLM family specified by the user via \code{fam}.
#'
#'   The function works by splitting the data into parts A and B,
#'   and computes a GLM fit on both parts.
#'   If \code{penalize == TRUE}, these fits use \code{cv.glmnet} from package \code{glmnet},
#'   otherwise they use \code{glmnet} with penalty set to 0. If
#'   \code{RP_function} (optional) is not supplied by
#'   the user, \code{randomForest} is used to predict remaining signal from the residuals from
#'   GLM fit on part A. The test statistic is proportional to the dot product
#'   between the random forest prediction and residuals from GLM fit on part B.
#'   If \code{nsplits} is greater than one, the above procedure is repeated \code{nsplits}
#'   times and the resulting p-values are aggregated using the approach from Meinshausen at al. (2012)
#'
#'   A user may supply their own residual prediction function to replace random forest via
#'   parameter \code{RP_function} (see Examples for use). The function must take as arguments
#'   an input matrix \code{XA}, vector \code{resA} (with length \code{nrow(XA)}) and matrix \code{XB}. Its
#'   role is to regress \code{resA} on input matrix \code{XA} with a preferred residual
#'   prediction method and output a vector with dimensions \code{nrow(XB)}
#'   that contains predictions of this fit on input \code{XB}.
#' @return If \code{output_all = FALSE}, the function outputs a single p-value.
#'   Otherwise it returns a list containing the aggregated p-value in \code{pval} and
#'   a vector of p-values from all splits in \code{pvals}.
#' @references
#' Janková, J., Shah, R. D., Bühlmann, P. and Samworth, R. (2019)
#' \emph{Goodness-of-fit testing in high-dimensional generalized linear models}
#' \url{https://arxiv.org/abs/1908.03606}
#' Meinshausen, N., Meier, L. and Bühlmann, P. (2012)
#' \emph{p-Values for High-Dimensional Regression}
#' Journal of the American Statistical Association, 104:488, 1671-1681
#' @seealso
#' @examples
#' # Testing for nonlinearity: Logistic link function
#'
#' set.seed(1)
#' X <- matrix(rnorm(300*30), 300, 30)
#' z <- X[, 1] + X[, 2]^4
#' pr <- 1/(1 + exp(-z))
#' y <- rbinom(nrow(X), 1, pr)
#' (out <- GRPtest(X, y, fam = "binomial", nsplits = 5))
#'
#' # Testing for nonlinearity: Define your own RP function
#' # use package xyz
#'
#' my_RP_function <- function(XA, resA, XB){
#'   xyz_fit <- xyz_regression(XA, resA)
#'   predict(xyz_fit, newdata = as.matrix(XB))[,5]
#' }
#'
#' library(xyz)
#' set.seed(2)
#' X <- matrix(rnorm(500*30), 500, 30)
#' z <- X[,1:3]%*%rep(1,3) + 1*X[, 1]*X[,5]
#' mu <- exp(z)
#' y <- rpois(n = nrow(X), lambda = mu)
#' (out <- GRPtest(X, y, fam = "poisson", RP_function = my_RP_function))
#'
#' @export
#' @import stats
#' @import glmnet
#' @import RPtests
#' @importFrom ranger ranger
#'
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
