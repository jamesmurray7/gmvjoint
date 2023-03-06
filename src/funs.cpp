#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;
static double const log2pi = std::log(2.0 * M_PI);

// Log-likelihoods --------------------------------------------------------
// Gaussian
double gaussian_ll(const arma::vec& Y,  const arma::vec& eta, const double sigma){ // important/misleading - sigma is the __variance__!!
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

// Gamma
//' @keywords internal
// [[Rcpp::export]]
double ll_Gamma(const arma::vec& Y, const double& shape, const arma::vec& mu){
  vec out = (shape - 1.) * log(Y) - lgamma(shape) - shape * log(mu) + 
    shape * log(shape) - (shape * Y)/mu;
  return sum(out);
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

//' @keywords internal
// [[Rcpp::export]]
double logfti(const arma::vec& b, const arma::rowvec& S, const arma::mat& SS, const arma::rowvec& Fi, const arma::mat& Fu,
              const double l0i, const arma::rowvec& haz, const int Delta, const arma::vec& gamma_rep, const arma::vec& zeta){
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  
  return as_scalar(
    temp + Delta * (S * zeta + Fi * (gamma_rep % b)) - haz * exp(SS * zeta + Fu * (gamma_rep % b))
  );
}

//' @keywords internal
// [[Rcpp::export]]
double joint_density(const arma::vec& b, const List Y, const List X, const List Z,                  // Longitudinal + Random effects.
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

  return -1.0 * (ll + as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + ll_Ti));;
}

// Quadrature - standard deviation of N(mu, tau^2).
//' @keywords internal
// [[Rcpp::export]]
List maketau(const List& S, const List& Z){
  int K = S.size();
  List out = List(K);
  for(int k = 0; k < K; k++){
    mat Zk = Z[k], Sk = S[k];
    out[k] = sqrt(diagvec(Zk * Sk * Zk.t()));
  }
  return out;
}

// Derivatives ------------------------------------------------------------
// Scores of the linear predictors eta = X(t)*beta+Z(t)*b;
arma::vec Score_eta_gauss(const arma::vec& eta, const arma::vec& Y, const double sigma){ // sigma = variance!
  int m = Y.size();
  mat V = mat(m, m, fill::eye);
  V *= sigma;
  return V.i() * (Y - eta);
}

// NB: Gaussian is __same__ with quadrature.

vec Score_eta_binom(const vec& eta, const vec& Y){
  return Y - exp(eta) / (exp(eta) + 1.0);
}

// Binomial d/d{eta} taken with quadrature.
vec Score_eta_binom_quad(const vec& eta, const vec& Y,
                         const vec& tau, const vec& w, const vec& v){
  int mi = Y.size(), gh = w.size();
  mat exp_part = mat(mi, gh);
  for(int l = 0; l < gh; l++){
    exp_part.col(l) = w[l] * exp(eta + tau * v[l]) / (exp(eta + tau * v[l]) + 1.);
  }
  return Y - sum(exp_part, 1);
}

vec Score_eta_poiss(const vec& eta, const vec& Y){
  return Y - exp(eta);
}

// Poisson d/d{eta} taken with quadrature
vec Score_eta_poiss_quad(const vec& eta, const vec& Y, const vec& tau){
  int mi = Y.size();
  vec Eexpmu = exp(eta + square(tau) * .5);
  return Y - Eexpmu;
}

vec Score_eta_genpois(const vec& eta, const vec& Y, const double phi, const mat& design){
  int q = design.n_cols;
  vec mu = exp(eta), grad = vec(q);
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(x + (Y - 1.) % (x % mu)/(mu + phi * Y) - x % mu / (phi + 1.));
  }
  return grad;
}

// GP1 d/d{eta} taken with quadrature
vec Score_eta_genpois_quad(const vec& eta, const vec& Y, const double phi, const mat& design,
                           const vec& tau, const vec& w, const vec& v){
  int q = design.n_cols, gh = w.size(), mi = Y.size();
  vec mu = exp(eta), grad = vec(q);
  vec exp_part1 = vec(mi), exp_part2 = vec(mi);
  // Work out parts we need to via quadrature.
  for(int l = 0; l < gh; l++){
    vec mu = exp(eta + tau * v[l]);
    exp_part1 += w[l] * mu / (phi * Y + mu);
    exp_part2 += w[l] * mu;
  }
  for(int qq = 0; qq < q; qq++){ // Work out column-by-column in {design}.
    vec x = design.col(qq);
    grad[qq] = sum(x + (Y - 1.) % (x % exp_part1) - x % exp_part2 / (phi + 1.));
  }
  return grad;
}

vec Score_eta_Gamma(const vec& eta, const vec& Y, const double shape, const mat& design){
  int q = design.n_cols;
  vec mu = exp(eta), grad = vec(q);
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(shape * x % Y/mu - shape * x);
  }
  return grad;
}

// TO DO:: Gamma quad version.

// Obtain kth derivative of log-likelihood wrt design matrix {X, Z}
// (Exploiting the structure of d/dbeta == d/db besides the design measure.)
vec get_long_score(const vec& eta, const vec& Y, const std::string family, const double sigma,
                   const mat& design){
  int p = design.n_cols;
  vec Score = vec(p);
  if(family == "poisson"){
    Score += design.t() * Score_eta_poiss(eta, Y);
  }else if(family == "gaussian"){
    Score += design.t() * Score_eta_gauss(eta, Y, sigma);
  }else if(family == "binomial"){
    Score += design.t() * Score_eta_binom(eta, Y);
  }else if(family == "genpois"){
    Score += Score_eta_genpois(eta, Y, sigma, design);
  }else if(family == "Gamma"){
    Score += Score_eta_Gamma(eta, Y, sigma, design);
  }
  return Score;
}

// For quadrature (Can't think of a neat way to do one function 25/11/22).
vec get_long_score_quad(const vec& eta, const vec& Y, const std::string family, const double sigma,
                        const mat& design, const vec& tau, const vec& w, const vec& v){
  int p = design.n_cols;
  vec Score = vec(p);
  if(family == "poisson"){
    Score += design.t() * Score_eta_poiss_quad(eta, Y, tau);
  }else if(family == "binomial"){
    Score += design.t() * Score_eta_binom_quad(eta, Y, tau, w, v);
  }else if(family == "genpois"){
    Score += Score_eta_genpois_quad(eta, Y, sigma, design, tau, w, v);
  }else if(family == "Gamma"){
    Score += Score_eta_Gamma(eta, Y, sigma, design);
  }
  return Score;
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec joint_density_ddb(const arma::vec& b, const List Y, const List X, const List Z,                  // Longitudinal + Random effects.
                            const arma::vec& beta, const arma::mat& D, const List sigma, const List family,
                            const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                            const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                            const arma::vec& gamma_rep, const arma::vec& zeta,
                            const List beta_inds, const List b_inds, const int K){
  vec Score = vec(b.size());
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    uvec b_k_inds = b_inds[k];
    vec beta_k = beta.elem(beta_k_inds); // Ensure indexing from zero!!
    vec b_k = b.elem(b_k_inds);
    vec eta = Xk * beta_k + Zk * b_k;
    Score.elem(b_k_inds) += get_long_score(eta, Yk, f, sigmak, Zk);
  }
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() % gamma_rep) - 
                          gamma_rep % (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma_rep % b)))));
}

// Second derivative of joint density is found by BFGS method in maximisation step.

//' @keywords internal
// [[Rcpp::export]]
arma::vec Sbeta(const arma::vec& beta, const List& X, const List& Y, const List& Z, const List& b, 
                const List& sigma, const List& family, const List& beta_inds, const int K,
                const bool quad, const List& tau, const arma::vec& w, const arma::vec& v){
  int p = beta.size();
  vec Score = vec(p);
  
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    uvec beta_k_inds = beta_inds[k];
    vec beta_k = beta.elem(beta_k_inds); 
    vec b_k = b[k];
    vec eta = Xk * beta_k + Zk * b_k;
    double sigmak = sigma[k];
//    if(f == "gaussian" || !quad || f != "binomial"){ // gaussian is the same, so do it here.
    if(f=="poisson"||f=="binomial"){
      vec tauk = tau[k]; // This is tau^2/2 (fron e.g. make tau2).
      Score.elem(beta_k_inds) += get_long_score_quad(eta, Yk, f, sigmak, Xk,
                 tauk, w, v);
      
    }else{
      Score.elem(beta_k_inds) += get_long_score(eta, Yk, f, sigmak, Xk);
    }
  }
  return Score;
}

// Hessian matrices for distns. of longit. part.
mat Hess_eta_gauss(const vec& eta, const vec& Y, const double sigma, const mat& design){
  int m = Y.size();
  mat V = mat(m, m, fill::eye);
  V *= sigma;
  return design.t() * -V.i() * design;
}

// Hess_eta_gauss is the __same__ with quadrature!  

mat Hess_eta_poiss(const vec& eta, const vec& Y, const mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec lambda = exp(eta);
  mat H = zeros<mat>(q, q);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (-1. * lambda[j]) * xj * xjT;
  }
  return H;
}

// Poisson d2/d{eta}2 taken with quadrature
arma::mat Hess_eta_poiss_quad(const arma::vec& eta, const arma::vec& Y, const arma::mat& design,
                              const arma::vec& tau){
  int mi = design.n_rows, q = design.n_cols;
  mat H = zeros<mat>(q, q);
  vec Eexpmu = exp(eta + square(tau) * .5);
  for(int j = 0; j < mi; j ++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H -= Eexpmu.at(j) * xj * xjT;
  }
  return H;
}

mat Hess_eta_binom(const vec& eta, const vec& Y, const mat& design){
  int mi = design.n_rows, q = design.n_cols;
  //vec mu = exp(eta);
  //vec denom = square(mu);
  vec expeta = exp(eta);
  vec cont = -1. * (expeta/(1. + expeta) - (expeta % expeta)/(square(1. + expeta)));
  mat H = zeros<mat>(q, q);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
//    H += (-mu[j]/(1. + denom[j])) * xj * xjT;
    H += cont.at(j) * xj * xjT;
  }
  return H;
}

// Binomial d2/d{eta}^2 taken with quadrature
mat Hess_eta_binom_quad(const arma::vec& eta, const arma::vec& Y, const arma::mat& design,
                        const arma::vec& tau, const arma::vec& w, const arma::vec& v){
  int mi = design.n_rows, q = design.n_cols, gh = w.size();
  mat H = zeros<mat>(q, q);
  vec exp_part = vec(mi), ones = vec(mi, fill::ones);
  for(int l = 0; l < gh; l++){
      vec exp_part_l = exp(eta + tau * v[l]);
      exp_part += w[l] * (exp_part_l/(ones + exp_part_l) % (ones - exp_part_l/(ones + exp_part_l)));
  }
  for(int j = 0; j < mi; j ++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H -= exp_part.at(j) * xj * xjT;
  }
  return H;
}

arma::mat Hess_eta_genpois(const arma::vec& eta, const arma::vec& Y, const double phi, const arma::mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec mu = exp(eta);
  mat H = zeros<mat>(q, q);
  vec inner_part = square(mu + phi * Y);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (phi * (Y[j] - 1.) * Y[j] / (inner_part[j]) - mu[j] / (phi + 1.)) * xj * xjT;
  }
  return H;
}

mat Hess_eta_genpois_quad(const vec& eta, const vec& Y, const double phi, const mat& design,
                          const vec& tau, const vec& w, const vec& v){
  int mi = design.n_rows, q = design.n_cols, gh = w.size();
  mat H = zeros<mat>(q, q);
  vec exp_part1a = vec(mi), exp_part1b = vec(mi), exp_part2 = vec(mi);
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + v[l] * tau;
    exp_part2 += w[l] * exp(eta_l);
    exp_part1a += w[l] * exp(eta_l);
    exp_part1b += w[l] * square(exp(eta_l) + phi * Y);
  }
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (phi * Y[j] * (Y[j] - 1.) * exp_part1a[j]/exp_part1b[j] - exp_part2[j] / (phi + 1.)) * xj * xjT;
  }
  return H;
}
                          
mat Hess_eta_Gamma(const vec& eta, const vec& Y, const double shape, const mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec mu = exp(eta);
  mat H = zeros<mat>(q, q);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (-shape * Y[j] / mu[j]) * xj * xjT;
  }
  return H;
}

// TODO:: Hessian for Gamma with quadrature.

mat get_long_hess(const vec& eta, const vec& Y, const std::string family, const double sigma, 
                  const mat& design){
  int p = design.n_cols;
  mat H = zeros<mat>(p, p);
  if(family == "poisson"){
    H += Hess_eta_poiss(eta, Y, design);
  }else if(family == "gaussian"){
    H += Hess_eta_gauss(eta, Y, sigma, design);
  }else if(family == "binomial"){
    H += Hess_eta_binom(eta, Y, design);
  }else if(family == "genpois"){
    H += Hess_eta_genpois(eta, Y, sigma, design);
  }else if(family == "Gamma"){
    H += Hess_eta_Gamma(eta, Y, sigma, design);
  }
  return H;
}

mat get_long_hess_quad(const vec& eta, const vec& Y, const std::string family, const double sigma,
                       const mat& design, const vec& tau, const vec& w, const vec& v){
  int p = design.n_cols;
  mat H = zeros<mat>(p, p);
  if(family == "poisson"){
    H += Hess_eta_poiss_quad(eta, Y, design, tau);
  }else if(family == "binomial"){
    H += Hess_eta_binom_quad(eta, Y, design, tau, w, v);
  }else if(family == "genpois"){
    H += Hess_eta_genpois_quad(eta, Y, sigma, design, tau, w, v);
  }else if(family == "Gamma"){
    H += Hess_eta_Gamma(eta, Y, sigma, design); //, tau2, w, v); // NYI
  }
  return H;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat Hbeta(const arma::vec& beta, const List& X, const List& Y, const List& Z, const List& b, 
                const List& sigma, const List& family, const List& beta_inds, const int K,
                const bool& quad, const List& tau, const arma::vec& w, const arma::vec& v){
  int P = beta.size();
  mat H = zeros<mat>(P, P);
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    int start = min(beta_k_inds), end = max(beta_k_inds);
    vec beta_k = beta.elem(beta_k_inds); // Ensure indexing from zero!!
    vec b_k = b[k];
    vec eta = Xk * beta_k + Zk * b_k;
//    if(f == "gaussian" || !quad || f != "binomial"){
    if(f == "poisson" || f == "binomial"){
      vec tauk = tau[k];
      H(span(start, end), span(start, end)) = get_long_hess_quad(eta, Yk, f, sigmak, Xk,
        tauk, w, v);
    }else{
      H(span(start, end), span(start, end)) = get_long_hess(eta, Yk, f, sigmak, Xk);
    }
  } // Return the (psuedo-) block diagonal.
  return H;
}

// Updates for dispersion parameters --------------------------------------

//' @keywords internal
// [[Rcpp::export]]
double vare_update(const arma::mat& X, const arma::vec& Y, const arma::mat& Z, const arma::vec& b, 
                   const arma::vec& beta, const arma::vec& tau, const arma::vec& w, const arma::vec& v){
  vec eta = X * beta + Z * b;
  int gh = w.size();
  double out = 0.0;
  for(int l = 0; l < gh; l++){
    vec rhs = Y - eta - tau * v[l];
    out += w[l] * as_scalar(rhs.t() * rhs);
  }
  return out;
}

//' @keywords internal
// [[Rcpp::export]]
List phi_update(const arma::vec& b, const arma::mat& X, const arma::vec& Y, const arma::mat& Z, 
                const arma::vec& beta, const double phi,
                const arma::vec& w, const arma::vec& v, const arma::vec& tau){
  int gh = w.size();
  vec eta = X * beta + Z * b;
  vec tau2 = tau;
  double rhs = sum(Y)/(1.+phi), lhs = sum(Y)/(pow(1.+phi,2.)), Score = 0., Hess = 0.;
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + tau2 * v[l];
    vec mu = exp(eta_l);
    Score += w[l] * sum(
      ((Y - 1.) % Y)/(mu + phi * Y) + (mu - Y)/(pow(phi + 1., 2.))
    );
    Hess += w[l] * sum(
      (2. * (Y - mu))/(pow(phi + 1., 3.)) - (square(Y) % (Y - 1.))/(square(mu + phi * Y))
    );
  }
  return List::create(_["Score"] = Score - rhs, _["Hessian"] = lhs + Hess);
}  

// Updates for the survival pair (gamma, zeta) ----------------------------
// Define the conditional expectation and then take Score AND Hessian via forward differencing
//' @keywords internal
// [[Rcpp::export]]
double Egammazeta(arma::vec& gammazeta, arma::vec& b, arma::mat& Sigma,
		   arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
		   arma::vec& haz, int Delta, arma::vec& w, arma::vec& v, Rcpp::List b_inds, int K){
	int gh = w.size(), q = Fu.n_cols; 
	vec g = gammazeta.head(K);
	vec z = gammazeta.subvec(K, gammazeta.size() - 1);
	vec gammas = vec(q);
	for(int k = 0; k < K; k++){
		double gk = g[k];
		uvec bk = b_inds[k];
		gammas.elem(bk) += gk;
	}
	rowvec gammas2 = gammas.t();
	mat Fugamma2 = Fu.each_row() % gammas2;
	mat A = Fugamma2 * Sigma * Fugamma2.t();
	vec mu = SS * z + Fu * (b % gammas);
	vec tau = sqrt(diagvec(A));
	double rhs = 0.;
	for(int l = 0; l < gh; l++){
		rhs += w[l] * as_scalar(haz.t() * exp(mu + v[l] * tau));
	}
	return as_scalar(Delta * (S * z + Fi * (b % gammas))) - rhs;
}


//' @keywords internal
// [[Rcpp::export]]
arma::vec Sgammazeta(arma::vec& gammazeta, arma::vec& b, arma::mat& Sigma,
                     arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                     arma::vec& haz, int Delta, arma::vec& w, arma::vec& v, Rcpp::List b_inds, int K, long double eps){
  int ps = gammazeta.size();
  vec out = vec(ps);
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
  for(int i = 0; i < ps; i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    double fdiff = Egammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K) - f0;
    out[i] = fdiff/(ge[i]-gammazeta[i]);
  }
  return out;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat Hgammazeta(arma::vec& gammazeta, arma::vec& b, arma::mat& Sigma,
                     arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                     arma::vec& haz, int Delta, arma::vec& w, arma::vec& v, Rcpp::List b_inds, int K, 
                     long double Seps, long double Heps){
  int ps = gammazeta.size();
  mat out = zeros<mat>(ps, ps);
  vec f0 = Sgammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, Seps);
  for(int i = 0; i < ps; i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * Heps;
    vec fdiff = Sgammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, Seps) - f0;
    out.col(i) = fdiff/(ge[i] - gammazeta[i]);
  }
  return 0.5 * (out + out.t());
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat lambdaUpdate(Rcpp::List survtimes, arma::mat& ft, arma::vec& gamma, arma::vec& zeta,
                       Rcpp::List S, Rcpp::List Sigma, Rcpp::List b, 
                       arma::vec& w, arma::vec& v, Rcpp::List b_inds){
  int gh = w.size(), n = S.size(), K = b_inds.size(), q = ft.n_cols;
  // Initialise store
  mat store = zeros<mat>(ft.n_rows, n);
  // Produce gamma vec
  vec gammas = vec(q);
  for(int k = 0; k < K; k++){   
    double gk = gamma.at(k);
    uvec bk = b_inds[k];
    gammas.elem(bk) += gk;
  }
  mat ftg = ft.each_row() % gammas.t();
  // Loop over subjects
  for(int i = 0; i < n; i++){
    vec survtimes_i = survtimes[i];
    mat Sigma_i = Sigma[i];
    vec b_i = b[i];
    rowvec S_i = S[i];
    mat A = ftg * Sigma_i * ftg.t();
    int ui = survtimes_i.size();
    // Loop over subject i's survived failure times.
    for(int j = 0; j < ui; j++){
      rowvec Fst = ft.row(j), Fstg = ftg.row(j);
      double mu = as_scalar(exp(S_i * zeta + Fst * (gammas % b_i)));
      double tau = sqrt(A.at(j,j));
      // Loop over GH nodes.
      for(int l = 0; l < gh; l++){
        store(j, i) += as_scalar(w[l] * mu * exp(v[l] * tau));
      }
    }
  }
  return store;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat joint_density_sdb(const arma::vec& b, const List Y, const List X, const List Z,                  // Longitudinal + Random effects.
                            const arma::vec& beta, const arma::mat& D, const List sigma, const List family,
                            const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                            const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                            const arma::vec& gamma_rep, const arma::vec& zeta,
                            const List beta_inds, const List b_inds, const int K, double eps){
  int n = b.size();
  mat out = zeros<mat>(n, n);
  vec f0 = joint_density_ddb(b, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
  for(int i = 0; i < n; i++){
    vec bb = b;
    double xi = std::max(b[i], 1.0);
    bb[i] = b[i] + (eps * xi);
    vec fdiff = joint_density_ddb(bb, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K) - f0;
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

/* *****
 * END-*
 * *****/
