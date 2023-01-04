#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
Omega.draw <- function(x){
  Omega.Var <- vcov(x)
  # Mean vector
  co <- x$coeffs
  sigma <- unlist(co$sig)
  sigma <- sigma[sigma != 0.0]
  Omega.mean <- setNames(
    c(
      c(vech(co$D)),
      co$beta, sigma,
      co$gamma, co$zeta
    ), names(x$SE))
  
  # Draw from N(Omega.mean, Omega.Var)
  draw <- MASS::mvrnorm(n = 1, mu = Omega.mean, Sigma = Omega.Var)
  
  # Re-construct Omega at this draw.
  D <- matrix(0, nrow(co$D), ncol(co$D)) # Need to check this works for K>1.
  D[lower.tri(D, T)] <- draw[grepl('^D\\[', names(draw))]
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  # Check this is pos-definite and transform if not.
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){ 
    D <- as.matrix(Matrix::nearPD(D)$mat)
  }
  beta <- draw[match(names(co$beta), names(draw))]
  gamma <- draw[match(names(co$gamma), names(draw))]
  zeta <- draw[match(names(co$zeta), names(draw))]
  .sigma <- draw[match(names(sigma), names(draw))]
  
  list(
    D = D, beta = beta, sigma = .sigma, gamma = gamma, zeta = zeta
  )
}
  