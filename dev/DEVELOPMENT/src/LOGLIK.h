#ifndef LOGLIK_H
#define LOGLIK_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

static double log2pi = log(2.0 * M_PI);

// Gaussian
double ll_Gaussian(const vec& mu, const vec& Y, const double& sigma){
  uword m_i = Y.size();
  vec S(m_i, fill::value(sqrt(sigma)));
  vec out = log_normpdf(Y, mu, S); // This is from Armadillo
  return sum(out);
}

// Poisson
double ll_Poisson(const vec& mu, const vec& Y){
  uword m_i = Y.size();
  vec out(m_i);
  for(uword i = 0; i < m_i; i++){
    out.at(i) = R::dpois(Y.at(i), mu.at(i), 1);
  }
  return sum(out);
}

// Binomial
double ll_Binomial(const vec& mu, const vec& Y){
  uword m_i = Y.size();
  vec out(m_i);
  for(uword i = 0; i < m_i; i++){
    out.at(i) = R::dbinom(Y.at(i), 1, mu.at(i), 1);
  }
  return sum(out);
}

// Gamma
double ll_Gamma(const vec& mu, const vec& Y, const vec& shape){
  vec scale = mu/shape;
  uword n = Y.size();
  vec out(n);
  for(uword i = 0; i < n; i++){
    out.at(i)  = R::dgamma(Y.at(i), shape.at(i), scale.at(i), 1);
  }
  return sum(out);
}

// Negative binomial
double ll_NegBin(const vec& mu, const vec& Y, const vec& phi){
  uword m_i = Y.size();
  vec out(m_i);
  for(uword i = 0; i < m_i; i++){
    out.at(i) = R::dnbinom_mu(Y.at(i), phi.at(i), mu.at(i), 1);
  }
  return sum(out);
}

// Generalised Poisson
// https://www.tandfonline.com/doi/pdf/10.1080/03610926.2011.564742?needAccess=true
double ll_GenPois(const vec& mu, const vec& Y, const vec& phi){
  vec frac = (mu + phi % Y) / (1. + phi);
  vec out = log(mu) + (Y - 1.) % log(mu + phi % Y) - lgamma(Y + 1.) - Y % log(1. + phi) - frac;
  return sum(out);
 }

// Multivariate Normal (for f(b|D))
// (see https://gallery.rcpp.org/articles/dmvnorm_arma/)
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

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


#endif 