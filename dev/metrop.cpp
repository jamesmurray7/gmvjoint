#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;
static double const log2pi = std::log(2.0 * M_PI);

double gaussian_ll(const arma::vec& Y,  const arma::vec& eta, const double sigma){ 
  int mi = Y.size();
  vec out = vec(mi);
  for(int i = 0; i < mi; i++){
    out[i] = R::dnorm(Y[i], eta[i], sqrt(sigma), 1);
  }
  return sum(out);
}

// Binomial
double binomial_ll(const arma::vec& Y, const arma::vec& eta){
  vec mu = exp(eta)/(1.0 + exp(eta));
  int mi = Y.size();
  vec out = vec(mi);
  for(int i = 0; i < mi; i++){
    out[i] = R::dbinom(Y[i], 1, mu[i], 1);
  }
  return sum(out);
}

// Poisson
double poisson_ll(const arma::vec& Y, const arma::vec& eta){
  int mi = Y.size();
  vec out = vec(mi);
  for(int i = 0; i < mi; i++){
    out[i] = R::dpois(Y[i], exp(eta[i]), 1);
  }
  return sum(out);
}


double logfti(const arma::vec& b, const arma::rowvec& S, const arma::mat& SS, const arma::rowvec& Fi, 
              const arma::mat& Fu, const double l0i, const arma::rowvec& haz, const int Delta, 
              const arma::vec& gamma_rep, const arma::vec& zeta){
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  
  return as_scalar(
    temp + Delta * (S * zeta + Fi * (gamma_rep % b)) - haz * exp(SS * zeta + Fu * (gamma_rep % b))
  );
}

double ll_Gamma2(const arma::vec& Y, const double& shape, const arma::vec& scale){
  int n = Y.size();
  vec out(n);
  for(int i = 0; i < n; i++){
    out[i]  = R::dgamma(Y[i], shape, scale[i], 1);
  }
  return sum(out);
}

// https://www.tandfonline.com/doi/pdf/10.1080/03610926.2011.564742?needAccess=true
//' @keywords internal
 // [[Rcpp::export]]
 double ll_genpois(const arma::vec& eta, const double phi, arma::vec& Y){
   vec mu = exp(eta);
   vec frac = (mu + phi * Y) / (1. + phi);
   vec out = log(mu) + (Y - 1.) % log(mu + phi * Y) - lgamma(Y + 1.) - Y * log(1. + phi) - frac;
   return sum(out);
 }



double joint_density(const arma::vec& b, const List Y, const List X, const List Z,
                     const arma::vec& beta, const arma::mat& D, const List sigma, const List family,
                     const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                     const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                     const arma::vec& gamma_rep, const arma::vec& zeta,
                     const List beta_inds, const List b_inds, const int K){
  double ll = 0.0; // Compile longitudinal parts ---------
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    uvec b_k_inds = b_inds[k];
    vec beta_k = beta.elem(beta_k_inds);  // ensure indices on R side are indexed from zero.
    vec b_k = b.elem(b_k_inds);
    vec eta = Xk * beta_k + Zk * b_k;
    if(f == "gaussian"){
      ll += gaussian_ll(Yk, eta, sigmak);
    }else if(f == "binomial"){
      ll += binomial_ll(Yk, eta);
    }else if(f == "poisson"){
      ll += poisson_ll(Yk, eta);
    }else if(f == "genpois"){
      ll += ll_genpois(eta, sigmak, Yk);
    }else if(f == "Gamma"){
      vec mu = exp(eta);
      ll += ll_Gamma2(Yk, sigmak, mu/sigmak); // == ll_Gamma(Yk, Sigmak, exp(eta)).
    }
  }
  // Rcout << "ll: " << ll << std::endl;
  int q = b.size();
  
  double ll_Ti = logfti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma_rep, zeta);
  
  return -1.0 * (ll + as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + ll_Ti));
}


// [[Rcpp::export]]
List metropolis(const arma::vec& b, const List Omega, const List Y, const List X, const List Z,
                const List family, const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                const arma::vec& gamma_rep, const List beta_inds, const List b_inds, const int K, const int q,
                const int burnin, const int N, const arma::mat& Sigma, const double tune){
  // Unpack Omega
  mat D = Omega["D"];
  List sigma = Omega["sigma"];
  vec beta = Omega["beta"];
  vec zeta = Omega["zeta"];
  // MC stuff
  int iters = burnin + N;
  mat out = mat(q, iters);
  out.col(0) += b;
  // Start
  int j = 1, num_accepts = 0;
  while(j < iters){
    double U = randu();
    vec b_current = out.col(j - 1);
    vec b_proposal = mvnrnd(b_current, tune * Sigma);
    double logf_current = -1. * joint_density(b_current, Y, X, Z, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double logf_proposal = -1. * joint_density(b_proposal, Y, X, Z, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double P = std::min(exp(logf_proposal - logf_current), 1.);
    if(U < P){
      b_current = b_proposal;
      if(j > burnin) num_accepts ++;
    }
    out.col(j) = b_current;
    j++ ;
  }
  
  // Remove burnin -----------
  out.shed_cols(0, burnin - 1);
  return List::create(_["walks"] = out,
                      _["burnin"] = burnin,
                      _["N"] = N,
                      _["AcceptanceRate"] = (double)num_accepts/(double)N);
}

double logdmvn(const vec& x, const vec& mu, const mat& S){
  vec xmmu = x - mu;
  return -0.5 * log2pi -0.5 * log(det(S)) - 0.5 * as_scalar(xmmu.t() * inv(S) * xmmu);
}

// [[Rcpp::export]]
List metropolis_hastings(const arma::vec& b, const List Omega, const List Y, const List X, const List Z,
                         const List family, const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                         const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                         const arma::vec& gamma_rep, const List beta_inds, const List b_inds, const int K, const int q,
                         const int burnin, const int N, const arma::mat& Sigma, const double tune){
  // Unpack Omega
  mat D = Omega["D"];
  List sigma = Omega["sigma"];
  vec beta = Omega["beta"];
  vec zeta = Omega["zeta"];
  // MC stuff
  int iters = burnin + N;
  mat out = mat(q, iters);
  out.col(0) += b;
  // Start
  int j = 1, num_accepts = 0;
  while(j < iters){
    double U = randu();
    vec b_current = out.col(j - 1);
    vec b_proposal = mvnrnd(b_current, tune * Sigma);
    double logf_current = -1. * joint_density(b_current, Y, X, Z, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double logf_proposal = -1. * joint_density(b_proposal, Y, X, Z, beta, D, sigma, family, Delta,
                                               S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double logQ_proposed = logdmvn(b_proposal, b_current, tune * Sigma);
    double logQ_current = logdmvn(b_current, b_proposal, tune * Sigma);
    double P = std::min(exp(logf_proposal + logQ_current - logf_current - logQ_proposed), 1.);
    if(U < P){
      b_current = b_proposal;
      if(j > burnin) num_accepts ++;
    }
    out.col(j) = b_current;
    j++ ;
  }
  
  // Remove burnin -----------
  out.shed_cols(0, burnin - 1);
  return List::create(_["walks"] = out,
                      _["burnin"] = burnin,
                      _["N"] = N,
                      _["AcceptanceRate"] = (double)num_accepts/(double)N);
}