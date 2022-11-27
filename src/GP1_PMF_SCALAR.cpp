#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Specifically for obtaining pmf in simulations (since requisite doesn't exist in R).
//' @keywords internal
// [[Rcpp::export]]
double GP1_pmf_scalar(const double mu, const double phi, const double Y){
  double L = mu * pow(mu + phi* Y, Y - 1.) * exp(-(mu+phi*Y)/(1.+phi)) / (pow(1. + phi, Y) * tgamma(Y + 1.));
  vec LL = vec(1, fill::value(L));
  LL.replace(datum::inf, 1e100);
  return as_scalar(LL);
}