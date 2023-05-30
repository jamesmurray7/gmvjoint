#include <RcppArmadillo.h>
#include <math.h>
#include "LOGLIK.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// vech(X) --> X ----------------------------------------------------------
arma::mat duplication_matrix(const int &n) {
  arma::mat out((n*(n+1))/2, n*n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      arma::vec u((n*(n+1))/2, arma::fill::zeros);
      u(j*n+i-((j+1)*j)/2) = 1.0;
      arma::mat T(n,n, arma::fill::zeros);
      T(i,j) = 1.0;
      T(j,i) = 1.0;
      out += u * arma::trans(arma::vectorise(T));
    }
  }
  return out.t();
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat vech2mat(const arma::vec& x, const int q){
  vec xx = duplication_matrix(q) * x;
  return reshape(xx, q, q);
}

// dmvn/dmvt --------------------------------------------------------------
//' @keywords internal
// [[Rcpp::export]]
double dmvn_fast(const arma::vec& x, const arma::vec& mean, const arma::mat& Sigma,
                 const bool log__ = true){
  return as_scalar(dmvnrm_arma_fast(x.t(), mean.t(), Sigma, log__));
}

//' @keywords internal
// [[Rcpp::export]]
double dmvt_fast(const arma::vec& x, const arma::vec& mean, const arma::mat& Sigma,
                 const double df, const bool log__ = true){
  return as_scalar(dmvt_arma_fast(x.t(), mean.t(), Sigma, df, log__));
}

// S(u)/S(t) --------------------------------------------------------------
//' @keywords internal
// [[Rcpp::export]]
double S_(const Rcpp::List& L, const arma::vec& gamma_rep, const arma::vec& zeta, const arma::vec& b){
  rowvec l0 = L["l0u"];
  mat SS = L["SS"];
  mat Fu = L["Fu"];
  return as_scalar(exp(-l0 * exp(SS * zeta + Fu * (gamma_rep % b))));
}