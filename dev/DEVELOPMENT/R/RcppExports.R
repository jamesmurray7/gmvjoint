# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
GP1_pmf_scalar <- function(mu, phi, Y) {
    .Call(`_DEVgmvjoint_GP1_pmf_scalar`, mu, phi, Y)
}

#' @keywords internal
vech2mat <- function(x, q) {
    .Call(`_DEVgmvjoint_vech2mat`, x, q)
}

#' @keywords internal
S_ <- function(L, gamma_rep, zeta, b) {
    .Call(`_DEVgmvjoint_S_`, L, gamma_rep, zeta, b)
}

make_eta <- function(X, Z, beta, b, beta_inds, b_inds) {
    .Call(`_DEVgmvjoint_make_eta`, X, Z, beta, b, beta_inds, b_inds)
}

make_tau <- function(Z, Sigma) {
    .Call(`_DEVgmvjoint_make_tau`, Z, Sigma)
}

#' @keywords internal
logfti <- function(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma_rep, zeta) {
    .Call(`_DEVgmvjoint_logfti`, b, S, SS, Fi, Fu, l0i, haz, Delta, gamma_rep, zeta)
}

joint_density <- function(b, Y, X, Z, W, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K) {
    .Call(`_DEVgmvjoint_joint_density`, b, Y, X, Z, W, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K)
}

joint_density_ddb <- function(b, Y, X, Z, W, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K) {
    .Call(`_DEVgmvjoint_joint_density_ddb`, b, Y, X, Z, W, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K)
}

Sbeta <- function(beta, X, Y, Z, W, b, sigma, family, beta_inds, b_inds, K, tau, w, v) {
    .Call(`_DEVgmvjoint_Sbeta`, beta, X, Y, Z, W, b, sigma, family, beta_inds, b_inds, K, tau, w, v)
}

Hbeta <- function(beta, X, Y, Z, W, b, sigma, family, beta_inds, b_inds, K, tau, w, v) {
    .Call(`_DEVgmvjoint_Hbeta`, beta, X, Y, Z, W, b, sigma, family, beta_inds, b_inds, K, tau, w, v)
}

appxE_Gammasigma <- function(sigma, eta, Y, tau, W, w, v) {
    .Call(`_DEVgmvjoint_appxE_Gammasigma`, sigma, eta, Y, tau, W, w, v)
}

appxE_NegBinsigma <- function(sigma, eta, Y, tau, W, w, v) {
    .Call(`_DEVgmvjoint_appxE_NegBinsigma`, sigma, eta, Y, tau, W, w, v)
}

appxE_GenPoissigma <- function(sigma, eta, Y, tau, W, w, v) {
    .Call(`_DEVgmvjoint_appxE_GenPoissigma`, sigma, eta, Y, tau, W, w, v)
}

sigma2_Gaussian_update <- function(eta, Y, tau, w, v) {
    .Call(`_DEVgmvjoint_sigma2_Gaussian_update`, eta, Y, tau, w, v)
}

#' @keywords internal, this assumes mu_surv, tau_surv not calculated prior.
lambda_hat <- function(b, Fu, SS, Sigma, gamma_rep, zeta, nev, w, v) {
    .Call(`_DEVgmvjoint_lambda_hat`, b, Fu, SS, Sigma, gamma_rep, zeta, nev, w, v)
}

#' @keywords internal
lambda_update <- function(b, Fu, SS, Sigma, survtimes, gamma_rep, zeta, nev, w, v) {
    .Call(`_DEVgmvjoint_lambda_update`, b, Fu, SS, Sigma, survtimes, gamma_rep, zeta, nev, w, v)
}

Egammazeta <- function(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, inds, K) {
    .Call(`_DEVgmvjoint_Egammazeta`, gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, inds, K)
}

#' @keywords internal
Sgammazeta <- function(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, eps) {
    .Call(`_DEVgmvjoint_Sgammazeta`, gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, eps)
}

#' @keywords internal
Hgammazeta <- function(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, eps) {
    .Call(`_DEVgmvjoint_Hgammazeta`, gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, eps)
}

#' @keywords internal
metropolis <- function(b, Omega, Y, X, Z, W, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, beta_inds, b_inds, K, q, burnin, N, Sigma, tune) {
    .Call(`_DEVgmvjoint_metropolis`, b, Omega, Y, X, Z, W, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, beta_inds, b_inds, K, q, burnin, N, Sigma, tune)
}

#' @keywords internal
dmvn_fast <- function(x, mn, Sigma, log__ = TRUE) {
    .Call(`_DEVgmvjoint_dmvn_fast`, x, mn, Sigma, log__)
}

#' @keywords internal
dmvt_fast <- function(x, mn, Sigma, df, log__ = TRUE) {
    .Call(`_DEVgmvjoint_dmvt_fast`, x, mn, Sigma, df, log__)
}

