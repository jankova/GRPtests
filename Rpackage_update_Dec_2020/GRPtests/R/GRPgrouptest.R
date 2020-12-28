#' Test significance of groups or individual predictors in high-dimensional generalized linear models
#'
#' The function can test significance of (potentially large) groups of predictors
#' in low- and high-dimensional generalized linear models.
#' Outputs a p-value.
#'
#' @param X Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param y Response vector.
#' @param fam Must be "gaussian", "binomial" or "poisson".
#' @param G A vector with indices of variables,
#'   whose significance we wish to ascertain, after controlling for variables in
#'   \code{X}. The size of \code{G} can be at most \code{p-2}.
#' @param B The number of bootstrap samples to approximate the distribution of
#'  the test statistic. Note that the p-value returned will always be
#'   at least \code{1/(B+1)}.
#' @param penalize If \code{TRUE}, penalization is used when fitting GLM models.
#' @details
#'   The function can test the significance of a set of variables in a generalized linear model,
#'   whose indices are specified by \code{G}.
#'   \code{penalize = TRUE} is needed for high-dimensional settings where the number of variables
#'   not in \code{G} is larger than the number of observations. We then employ a penalized regression
#'   to regress \code{y} on to these variables implemented in \code{cv.glmnet} from package \code{glmnet}.
#'   For the low-dimensional case, an unpenalized regression may be used.
#' @return The output is a single p-value.
#' @references Janková, J., Shah, R. D., Bühlmann, P. and Samworth, R. (2019)
#' \emph{Goodness-of-fit testing in high-dimensional generalized linear models}
#' \url{https://arxiv.org/abs/1908.03606}
#' @seealso
#' @examples
#' # Testing significance of a group of predictors in logistic regression
#'
#' set.seed(1)
#' X <- matrix(rnorm(300*50), 300, 50)
#' z <- X[, 1:5] %*% rep(1, 5)
#' pr <- 1/(1 + exp(-z))
#' y <- rbinom(nrow(X), 1, pr)
#' (out <- GRPgrouptest(X, y, fam = "binomial", G = 5:10, B = 1000))
#'
#' @export
#' @import stats
#' @import randomForest
#' @import glmnet
#' @import RPtests
#'
#'

GRPgrouptest <- function(X, y, fam = c("gaussian", "binomial", "poisson"),
                         G, B = 1000L, penalize = ifelse(p-length(G) >= n, TRUE, FALSE)){

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


  fam <- match.arg(fam)

  # if(fam != "gaussian" && fam != "binomial" && fam != "poisson"){
  #   stop("family must be gaussian, binomial or poisson; try changing parameter fam.")
  # }

  if (!is.numeric(G)){
    stop("G must be a numeric vector.")
  } else if (min(G) <= 0 || max(G) > p){
    stop("G must contain integers from the set {1, 2, ..., p}.")
  } else if (sum(duplicated(G)) > 0){
    stop("No duplicate values are allowed in G.")
  } else if (length(G) > p - 2){
    stop("The size of G can be at most p-2.")
  } else if (sum(round(G) != G) > 0){
    stop("G must be a vector of integers.")
  }

  if ((p - length(G) >= n-1) && penalize == FALSE) {
    stop("When penalize=FALSE we must have ncol(X) < nrow(X)-1; try setting penalize=TRUE.")
  }

  if (!is.numeric(B) || round(B) != B){
    stop("B must be an integer.")
  }

  pval_computed <- group_test(X, y, G, fam, B, penalize)

  return(pval_computed)
}
