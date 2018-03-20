
#' cv_cdmm
#'
#' Chooses lambda for the compositional lasso by cross-validated prediction error
#'
#' @param y Vector of outcomes
#' @param x Matrix of log relative abundances
#' @param nfolds Number of folds for cross validation
#' @param lamseq Sequence of lambda values to try
#' @param nlam Number of lambda values to try; default 100
#' @param rlam Minimum lambda value as a fraction of maxlam; default 1/nlam
#' @param mu Lagrange scale parameter; default 1
#' @param std Logical indicating whether data should be standardized; default TRUE
#' @param maxv Maximum number of nonzero variables; default = 0.4*n
#' @param maxit Maximum number of iterations (two measures of convergence); default = c(20, 50)
#' @param tol Tolerance for convergence (same two measures); default = c(1e-4, 1e-7)
#' @export
#'
#'
cv_cdmm <- function(y, x, nfolds, lamseq, nlam=100, rlam=1/nlam, mu=1, std=TRUE, maxv=0.4*length(y), maxit=c(20,50), tol=c(1e-4, 1e-7)) {
  reord <- sample(1:nrow(x))
  folds <- sapply(1:nfolds, FUN = function(i) reord[which( (1:nrow(x) %% nfolds + 1) == i)])
  if (missing(lamseq)) {
      cl.res <- cdmm(y = y, x = x, nlam = nlam, rlam, mu, std, maxv, maxit, tol)  ## initial fit
      lamseq <- cl.res$lamseq
  } else {
      nlam = length(lamseq)
      cl.res <- cdmm(y = y, x = x, lamseq = lamseq, mu, std, maxv, maxit, tol)  ## initial fit
  }

  cl.pe <- matrix(nrow = nrow(x), ncol = nlam)   ## nlam = 100

  for (i in 1:nfolds) {
    x.train <- x[-folds[[i]], ]
    x.test <- x[folds[[i]], ]
    y.train <- y[-folds[[i]]]
    y.test <- y[folds[[i]]]

    this.res <- cdmm(y = y.train, x = x.train, lamseq = lamseq, mu, std, maxv, maxit, tol)
    intmat <- matrix(rep(this.res$int, nrow(x.test)), nrow = nrow(x.test), byrow = TRUE)
    cl.pe[folds[[i]], ] <- y.test - x.test %*% this.res$sol - intmat
  }

  cl.pe.res <- apply(cl.pe, 2, FUN = function(xx) sum(xx^2/nrow(x)))
  #plot(cl.pe.res ~ lamseq)
  cl.beta <- cl.res$sol[, which.min(cl.pe.res)]
  cl.int <- cl.res$int[which.min(cl.pe.res)]
  return(list(cv.beta = cl.beta, cv.int = cl.int, pe = cl.pe, fullfit = cl.res))
}



