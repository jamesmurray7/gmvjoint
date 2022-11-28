long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ bin
family <- list('genpois')
data <- simData(beta = t(c(2, -0.1, 0.1, -0.2)), sigma = c(1.),
                family = family, gamma = .5,
                D = diag(c(.25, .09)))
data <- data$data
control <- list()

## run up to start EM in joint.R
## And to after sigma.update in EMupdate

Ellgp <- function(phi, y, X, z, b, tau, tau2){
  eta <- X[[1]]%*%beta[1:4]+z[[1]]%*%b[1:2]
  y <- y[[1]]
  gh <- length(w); mi <- length(y)
  p1 <- p2 <- matrix(0, nr = mi, nc = gh)
  for(l in 1:gh){
    p1[,l] <- w[l] * log(exp(eta + tau[[1]] * v[l]) + phi * y)
    p2[,l] <- w[l] * exp(eta + tau2[[1]] * v[l])
  }
  p1 <- rowSums(p1); p2 <- rowSums(p2)
  print(p1)
  print(p2)
  sum(
    (y-1) * p1 - y * log(1+phi) - (p2 + phi * y)/(1+phi)
  )
}

l <- S <- H <- numeric(n)
for(i in 1:n){
  l[i] <- Ellgp(sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                tau = tau[[i]], tau2 = tau2[[i]])
  S[i] <-  pracma::grad(Ellgp, sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                        tau = tau[[i]], tau2 = tau2[[i]])
  H[i] <-  pracma::hessian(Ellgp, sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                           tau = tau[[i]], tau2 = tau2[[i]])
}

Ellgp2 <- function(phi, y, X, z, b, tau, tau2, S){
  eta <- X[[1]]%*%beta[1:4]+z[[1]]%*%b[1:2]
  y <- y[[1]]
  gh <- length(w); mi <- length(y)
  p1 <- p2 <- numeric(mi)
  for(j in 1:mi){
    p1[j] <- integrate(function(x) log(exp(x) + sigma[[1]] * y[j])*dnorm(x, mean = eta[j], sd = tau[j]),
                       lower = -qnorm(.99, eta[j], tau[j]), upper = qnorm(.99, eta[j], tau[j]))$value # hardcode range
    p2[j] <- integrate(function(x) exp(x) * dnorm(x, mean = eta[j], sd = tau[j]),
                       lower = -qnorm(.99, eta[j], tau[j]), upper = qnorm(.99, eta[j], tau[j]))$value # hardcode range
  }
  print(p1)
  print(p2)
  sum(
    (y-1) * p1 - y * log(1+phi) - (p2 + phi * y)/(1+phi)
  )
}

l2 <- S2 <- H2 <- numeric(n)
for(i in 1:n){
  l2[i] <- Ellgp2(sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                tau = tau[[i]], tau2 = tau2[[i]], S = Sigma[[i]])
  S2[i] <-  pracma::grad(Ellgp2, sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                        tau = tau[[i]], tau2 = tau2[[i]], S = Sigma[[i]])
  H2[i] <-  pracma::hessian(Ellgp2, sigma[[1]]  , y = Y[[i]], X = X[[i]], z = Z[[i]], b = b.hat[[i]],
                           tau = tau[[i]], tau2 = tau2[[i]], S = Sigma[[i]])
}

plot(l, l2);abline(0,1)

sigma[[1]] - sum(S)/sum(H)
sigma[[1]] - sum(S2)/sum(H2)
