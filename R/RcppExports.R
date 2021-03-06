# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' cdmm_amp
#'
#' Runs Lin's compositional lasso
#'
#' @param y Outcome vector
#' @param x OTU matrix
#' @param fac Vector of scaling factors for OTU matrix
#' @param lams Sequence of lambda values to try
#' @param mu Lagrange scale parameter
#' @param maxv Max number of nonzero variables
#' @param maxit Max number of iterations allowed
#' @param tol Parameter defining convergence
#' @export
cdmm_amp <- function(y, x, fac, lams, mu, maxv, maxit, tol) {
    .Call(`_CompLasso_cdmm_amp`, y, x, fac, lams, mu, maxv, maxit, tol)
}

#' signC
#'
#' Finds the sign of a double
#'
#' @param x Number to find the sign of
#' @export
signC <- function(x) {
    .Call(`_CompLasso_signC`, x)
}

#' absC
#'
#' Finds the absolute value of a double
#'
#' @param x Number to take abs value of
#' @export
absC <- function(x) {
    .Call(`_CompLasso_absC`, x)
}

#' absC2
#'
#' Finds the absolute value of each element in a numeric vector
#'
#' @param x Vector to take abs value of
#' @export
absC2 <- function(x) {
    .Call(`_CompLasso_absC2`, x)
}

#' softC
#'
#' Thresholds x by lambda
#'
#' @param x Number to soft threshold
#' @param lam Parameter by which to soft threshold
#' @export
softC <- function(x, lam) {
    .Call(`_CompLasso_softC`, x, lam)
}

#' softC2
#'
#' Soft thresholds each element of a vector x by lambda
#'
#' @param x Vector whose elements should be soft thresholded
#' @param lam Parameter by which to soft threshold
#' @export
softC2 <- function(x, lam) {
    .Call(`_CompLasso_softC2`, x, lam)
}

#' sqrtC
#'
#' Finds the square root of a double
#'
#' @param x Number to take sqrt of
#' @export
sqrtC <- function(x) {
    .Call(`_CompLasso_sqrtC`, x)
}

