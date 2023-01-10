#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double log2pi = log(2.0 * M_PI); 

// Fast multivariate normal and t density calculations (see https://gallery.rcpp.org/articles/dmvnorm_arma/)
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = true) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dmvt_arma_fast(arma::mat const &x,             // fast DMVT inspired by the above.
                         arma::rowvec const &mean,
                         arma::mat const &sigma, 
                         double const df,
                         bool const logd = true){
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = lgamma(((double)xdim + df)/2.0) - lgamma(df/2.0) - (double)xdim/2.0 * log(M_PI * df), 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * (df + (double)xdim) * log(arma::dot(z, z)/df + 1.0);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

//' @keywords internal
// [[Rcpp::export]]
double S_(const Rcpp::List& L, const arma::vec& gamma_rep, const arma::vec& zeta, const arma::vec& b){
	rowvec l0 = L["l0u"];
	mat SS = L["SS"];
	mat Fu = L["Fu"];
	return as_scalar(exp(-l0 * exp(SS * zeta + Fu * (gamma_rep % b))));
}
