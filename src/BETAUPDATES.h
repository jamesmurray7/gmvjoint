#ifndef BETAUPDATES_H
#define BETAUPDATES_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace arma;

// Scores with respect to eta ---------------------------------------------
// d/deta with 

// Gaussian
vec d_eta_Gaussian(const vec& eta, const vec& Y, const double& sigma){
  uword m_i = Y.size();
  mat V(m_i, m_i, fill::eye);
  V *= sigma;
  return V.i() * (Y - eta);
}

mat d2_eta_Gaussian(const vec& Y, const double& sigma){
  uword m_i = Y.size();
  mat V(m_i, m_i, fill::eye);
  V *= sigma;
  return -V.i();
}

// Poisson -- Log-normal, evaluate at median (since distn tends to be skew.)
vec d_eta_Poisson(const vec& eta, const vec& Y){
  return Y - trunc_exp(eta);
}

vec d2_eta_Poisson(const vec& eta, const vec& Y){
  return -1. * trunc_exp(eta);
}

// Binomial
vec d_eta_Binomial(const vec& eta, const vec& Y, const vec& tau,
                   const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp(eta.size());
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * trunc_exp(this_eta)/(1. + trunc_exp(this_eta));
  }
  return Y - Exp;
}

vec d2_eta_Binomial(const vec& eta, const vec& tau, const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp(eta.size());
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * trunc_exp(this_eta) % pow(trunc_exp(this_eta) + 1., -2.0);
  }
  return -Exp;
}

// Gamma
vec d_eta_Gamma(const vec& eta, const vec& Y, const vec& shape, const vec& tau,
                const vec& w, const vec& v){
  uword gh = w.size(), n = eta.size();
  vec out(n);
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    out += w[l] * (shape % (Y - trunc_exp(this_eta)) / trunc_exp(this_eta));
  }
  return out;
}

vec d2_eta_Gamma(const vec& eta, const vec& Y, const vec& shape, const vec& tau,
                 const vec& w, const vec& v){
  uword gh = w.size(), n = eta.size();
  vec out(n);
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    out += w[l] * shape % (-1. * trunc_exp(-1. * this_eta)) % Y;
  }
  return out;
}

// Negative binomial
vec d_eta_NegBin(const vec& eta, const vec& Y, const vec& phi, const vec& tau,
                 const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp(eta.size());
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * trunc_exp(this_eta)/(phi + trunc_exp(this_eta));
  }
  return Y - (Y + phi) % Exp;
}

vec d2_eta_NegBin(const vec& eta, const vec& Y, const vec& phi, const vec& tau,
                  const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp(eta.size());
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * phi % trunc_exp(this_eta) % pow(phi + trunc_exp(this_eta), -2.);
  }
  return -1. * (Y + phi) % Exp;
}

// Generalised Poisson
vec d_eta_GenPois(const vec& eta, const vec& Y, const vec& phi, const vec& tau,
                  const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp1(eta.size()), Exp2 = trunc_exp(eta + square(tau)/2.);
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau* v[l];
    Exp1 += w[l] * trunc_exp(this_eta)/(trunc_exp(this_eta) + phi % Y);
  }
  return 1. + (Y - 1.) % Exp1 - Exp2;
}

vec d2_eta_GenPois(const vec& eta, const vec& Y, const vec& phi, const vec& tau,
                   const vec& w, const vec& v){
  uword gh = w.size();
  vec Exp1(eta.size()), Exp2 = trunc_exp(eta + square(tau)/2.);
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau* v[l];
    Exp1 += w[l] * phi % Y % trunc_exp(this_eta) % pow(trunc_exp(this_eta) + phi % Y, -2.);
  }
  return (Y - 1.) % Exp1 - Exp2;
}

mat form_hess(const vec& d2, const mat& X){
  mat dm = diagmat(d2);
  return X.t() * dm * X;
}

#endif