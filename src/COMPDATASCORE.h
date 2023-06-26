#ifndef COMPDATASCORE_H
#define COMPDATASCORE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

// Scores with respect to eta ---------------------------------------------
// Some score/hessians for beta is then just X.t() * Score_...
// Gradient for b in optim step is always Z.t() * Score_...

// Gaussian
vec score_Gaussian(const vec& eta, const vec& Y, const double& sigma){
  uword m_i = Y.size();
  mat V(m_i, m_i, fill::eye);
  V *= sigma;
  return V.i() * (Y - eta);
}

// Poisson
vec score_Poisson(const vec& eta, const vec& Y){
  return Y - trunc_exp(eta);
}

// Binomial
vec score_Binomial(const vec& eta, const vec& Y){
  return Y - trunc_exp(eta)/(1. + trunc_exp(eta));
}

// Gamma
vec score_Gamma(const vec& eta, const vec& Y, const vec& shape){
  // return Y % shape % trunc_exp(-eta) - shape;
  vec neg_mu = trunc_exp(-1. * eta);
  // return mu % (Y/mu - 1.) % shape/mu;
  return shape % (neg_mu % Y - 1.);
}

// Negative binomial
vec score_NegBin(const vec& eta, const vec& Y, const vec& phi){
  vec mu = trunc_exp(eta);
  return Y - mu % (Y + phi) / (mu + phi);
}

// Generalised Poisson
vec score_GenPois(const vec& eta, const vec& Y, const vec& phi){
  vec mu = trunc_exp(eta);
  return 1. + mu % (Y - 1)/(mu + phi % Y) - mu/(1. + phi);
}


// Multivariate Normal (for f(b|D))
vec grad_b(const vec& b, const mat& D){
  return -D.i() * b;
}

// Derivative of log f(T_i|b, ...) wrt b
vec Score_Ti(const arma::vec& b, const arma::rowvec& S, const arma::mat& SS, const arma::rowvec& Fi, const arma::mat& Fu,
             const double l0i, const arma::rowvec& haz, const int Delta, const arma::vec& gamma_rep, const arma::vec& zeta){
  vec expo = exp(SS * zeta + Fu * (gamma_rep % b));
  return Delta * (Fi.t() % gamma_rep) - gamma_rep % (Fu.t() * (haz.t() % expo));
}
  

#endif 