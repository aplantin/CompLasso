#' cdmm
#'
#' Wrapper function for C++ function that runs the compositional lasso
#' @useDynLib CompLasso
#' @importFrom Rcpp sourceCpp
#'
#' @param y Vector of outcomes
#' @param x Matrix of log relative abundances
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
cdmm <- function(y, x, lamseq, nlam=100, rlam=1/nlam, mu=1, std=TRUE, maxv=0.4*length(y), maxit=c(20, 50), tol=c(1e-4, 1e-7)) {
  if (std == TRUE) {
    y <- scale(y, scale=FALSE)
    x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
    fac <- 1/attr(x, "scaled:scale")
  } else {
    fac <- rep(1, ncol(x))
  }
  if (missing(lamseq)) {
    lam.max <- max(abs(crossprod(x, y)))
    lamseq <- lam.max*exp(seq(0, log(rlam), length=nlam))
  }
  sol1 <- cdmm_amp(y = y, x = x, fac = fac, lams = lamseq,
                   mu = mu, maxv = maxv, as.integer(maxit), tol = tol)
  sol <- matrix(sol1, nrow = ncol(x), ncol = nlam, byrow = FALSE)
  sol <- sol*fac
  if (std == TRUE) {
    int <- attr(y, "scaled:center") - drop(crossprod(sol, attr(x, "scaled:center")))
  }
  return(list(sol = sol, int = int, lamseq = lamseq))
}


#' gic.cdmm
#'
#' Wrapper function for C++ function that runs the compositional lasso,
#' using the generalized information criterion (GIC) to select lambda
#'
#' @param y Vector of outcomes
#' @param x Matrix of log relative abundances
#' @param nlam Number of lambda values to try; default 100
#' @export
#'
gic.cdmm <- function(y, x, nlam) {
  res <- cdmm(y, x, nlam)
  n <- length(y)
  p <- ncol(x)
  fit <- matrix(res$int, nrow = n, ncol = nlam, byrow=TRUE) + x %*% res$sol
  lsig2.lam <- log(apply(fit, 2, FUN = function(ff) sum((y - ff)^2)/n))
  ft <- log(log(n))*log(max(p, n))/n
  s.lam <- colSums(res$sol != 0)
  gic <- lsig2.lam + ft*(s.lam - 1)
  ilam <- which.min(gic)
  return( list(bet = res$sol[, ilam], lam = res$lam[ilam], int = res$int[ilam] ))
}

