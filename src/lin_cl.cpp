#include "service.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' cdmm_amp
//'
//' Runs Lin's compositional lasso
//'
//' @param y Outcome vector
//' @param x OTU matrix
//' @param fac Vector of scaling factors for OTU matrix
//' @param lams Sequence of lambda values to try
//' @param mu Lagrange scale parameter
//' @param maxv Max number of nonzero variables
//' @param maxit Max number of iterations allowed
//' @param tol Parameter defining convergence
//' @export
// [[Rcpp::export]]
SEXP cdmm_amp(NumericVector y, NumericMatrix x, NumericVector fac, NumericVector lams,
	      double mu, double maxv, IntegerVector maxit, NumericVector tol) {
  double n = y.length();    // i
  double p = x.ncol();      // j
  int nlam = lams.length(); // k
  int maxit_cd = maxit[0];  // itcd
  int maxit_mm = maxit[1];  // itmm
  double tol_cd = tol[0];
  double tol_mm = tol[1];

  NumericVector sol(p*nlam);
  NumericVector bet(p);
  NumericVector res(n);

  for (int i = 0; i < p; i++) bet[i] = 0.0;
  for (int i = 0; i < n; i++) res[i] = y[i];
  double mu1 = 1.0/(1.0 + mu); double sum = 0.0; double norm = 0.0;
  for (int j = 0; j < nlam; j++) {
    double thislam = lams[j]; double lam1 = thislam*mu1; double alp = 0.0;
    for (int itmm = 0; itmm < maxit_mm; itmm++) {
      for (int itcd = 0; itcd < maxit_cd; itcd++) {
        double del = 0.0; double norm_old = norm; norm = 0.0;
        for (int i = 0; i < p; i++) {
          double old = bet[i]; double t0 = 0.0;
          for (int k = 0; k < n; k++) t0 += res[k]*x[k + i*n];
          t0 = (t0 - mu*(alp + sum))*mu1 + bet[i];
          bet[i] = softC(t0, lam1);
          if (bet[i] != old) {
            double dif = bet[i] - old; sum += fac[i]*dif;
            for (int k = 0; k < n; k++) res[k] -= x[k + i*n]*dif;
            del += std::abs(dif);
          }
          norm += std::abs(bet[i]);
        }
        if (del < tol_cd*norm_old) break;
      }
      if (std::abs(sum) < tol_mm) break;
      alp += sum;
    }
    int nvar = 0;
    for (int i = 0; i < p; i++) {
      sol[i + j*p] = bet[i];
      if (bet[i] != 0.0) nvar++;
    }
    if (R_FINITE(maxv) && nvar > maxv) break;
  }

  //  return List::create(Named("sol") = sol, Named("lams") = lams);
  return sol;
}

