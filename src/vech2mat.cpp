#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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