#' @keywords internal
E_shape.b <- function(shape, X, Y, Z, tau, beta, b, w, v){ # Shape for gamma via pracma::grad/hessian
  shape <- shape
  eta <- X %*% beta + Z %*% b
  out <- numeric(length(w))
  for(l in 1:length(w)){
    mu.l <- exp(eta + tau * v[l])
    out[l] <- w[l] * ll_Gamma(Y, shape, mu.l)
  }
  sum(out)
}

#' @keywords internal
#' @importFrom pracma grad hessian
shape_update <- function(shape, X, Y, Z, tau, beta, b, w, v){
  out <- setNames(vector('list', 2), c("Score", "Hessian"))
  out$Score <- pracma::grad(E_shape.b, shape, X = X, Y = Y, Z = Z, tau = tau, beta = beta, b = b, w = w, v = v)
  out$Hessian <- pracma::hessian(E_shape.b, shape, X = X, Y = Y, Z = Z, tau = tau, beta = beta, b = b, w = w, v = v)
  out
}