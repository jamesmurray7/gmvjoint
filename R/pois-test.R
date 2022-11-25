beta <- do.call(rbind, replicate(1, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
gamma <- c(0.3)
D <- diag(c(0.25, 0.09))
data <- simData(ntms = 15, beta = beta, D = D, 
                family = list('poisson'), zeta = c(0, -0.2),
                sigma = c(0.16), gamma = gamma)$data

# Specify formulae and target families
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ bin
family <- list('poisson')

fit <- joint(long.formulas, surv.formula, data, family)


tau2 <- mapply(function(Z, S) diag(tcrossprod(Z[[1]] %*% S, Z[[1]]))/2, Z = Z, S = Sigma)

Sb <- mapply(function(X, Y, Z, b, tau2){
  mu <- X[[1]] %*% beta + Z[[1]] %*% b
  expbit <- matrix(0, nr = length(mu), nc = 3)
  for(l in 1:3) expbit[,l] <- w[l] * exp(mu + tau2*v[l])
  -crossprod(X[[1]], rowSums(expbit)) + c(crossprod(Y[[1]], X[[1]]))
}, X = X, Y = Y, Z = Z, b = b.hat, tau2 = tau2,SIMPLIFY = F)

Hb <- mapply(function(X, Y, Z, b, tau2){
  mu <- X[[1]] %*% beta + Z[[1]] %*% b
  expbit <- matrix(0, nr = length(mu), nc = 3)
  for(l in 1:3) expbit[,l] <- w[l] * exp(mu + tau2*v[l])
  expbit <- rowSums(expbit)
  H <- matrix(0, length(beta), length(beta))
  for(j in 1:length(mu)){
    Xj <- X[[1]][j,,drop=T]
    H <- H + expbit[j] * (tcrossprod(Xj))
  }
  -H
}, X = X, Y = Y, Z = Z, b = b.hat, tau2 = tau2,SIMPLIFY = F)