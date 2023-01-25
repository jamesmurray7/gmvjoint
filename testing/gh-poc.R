# Simulate some data ------------------------------------------------------
# say Poisson
rm(list=ls())
pois <- simData(family = list('poisson'), sigma = .1, beta = t(c(1,-0.1,0.1,0.2)),
                gamma = 0.5, zeta = c(0, -0.15))
data <- pois$data
long.formulas <- list(Y.1 ~ time + cont + bin + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ bin
family <- list('poisson')


# Parse a la joint --------------------------------------------------------
# Initial parsing ----
formulas <- lapply(long.formulas, parseFormula)
surv <- parseCoxph(surv.formula, data)
n <- surv$n; K <- length(family)

# Initial conditons ----
inits.long <- Longit.inits(long.formulas, data, family)
# Suss out indices of b_k and beta_k.
b.inds <- lapply(1:K, function(k){
  nm <- inits.long$responses[k]
  which(grepl(nm, colnames(inits.long$b)))
})
beta.inds <- lapply(1:K, function(k){
  nm <- inits.long$responses[k]
  which(grepl(nm, names(inits.long$beta.init)))
})
q <- length(do.call(c, b.inds))

inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas, b.inds)

# Longitudinal parameters
beta <- inits.long$beta.init
D <- inits.long$D.init
sigma <- inits.long$sigma.init # dispersion / resid. variance / 0 otherwise.
b <- lapply(1:n, function(i) inits.long$b[i, ])
# Survival parameters
zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
names(zeta) <- paste0('zeta_', names(zeta))
gamma <- inits.surv$inits[grepl('gamma\\_', names(inits.surv$inits))]

# Longitudinal and survival data objects ----
dmats <- createDataMatrices(data, formulas)
sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)

X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
m <- lapply(Y, function(y) sapply(y, length))
# survival
Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
l0 <- sv$l0
S <- sv$S; SS <- sv$SS

# Assign family to joint density and parameter updates ----
# Ensure family is a string, not a function i.e. if user supplies `Gamma` instead of `"Gamma"`.
family <- lapply(family, function(family) if("function"%in%class(family)) family()$family else family)

Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
            beta, unlist(sigma)[inits.long$sigma.include], gamma, zeta)
# Temporary E-step --------------------------------------------------------

# Unpack Omega, the parameter vector
D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) # Indexed for C++ use.
gamma.rep <- rep(gamma, sapply(b.inds, length))

# b.hat
b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
  optim(b, joint_density, joint_density_ddb,
        Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
        Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
        beta_inds = beta.inds2, b_inds = b.inds2, K = K,
        method = 'BFGS', hessian = T)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
Fu = Fu, l0u = l0u, SIMPLIFY = F)

b.hat <- lapply(b.update, el, 1)
Sigma <- lapply(b.update, function(x) solve(x$hess))


# gauss.quad --------------------------------------------------------------
GH1 <- statmod::gauss.quad(3, 'hermite')
w <- GH1$we; v <- GH1$no

Eb1 <- mapply(function(mu, S){ # This should simply === b.hat!
  sig <- sqrt(diag(S))
  out <- matrix(NA, nr = 3, nc = 2)
  for(l in 1:3)
    out[l,] <- w[l] * (sqrt(2) * sig * v[l] + mu)
  colSums(out)/sqrt(pi)
}, mu = b.hat, S = Sigma, SIMPLIFY = F)

eta <- mapply(function(X, Z, b) X[[1]] %*% beta + Z[[1]] %*% b, X = X, Z = Z, b = b.hat, SIMPLIFY = F)
tau <- mapply(function(Z, S) sqrt(2) * sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]]))), Z = Z, S = Sigma, SIMPLIFY = F)
Eexpeta1 <- mapply(function(eta, tau){
  out <- matrix(NA, nr = length(eta), nc = 3)
  for(l in 1:3)
    out[, l] <- w[l] * exp(tau * v[l] + eta)
  rowSums(out)/sqrt(pi)
}, eta = eta, tau = tau, SIMPLIFY = F)


# gauss.quad.prob ---------------------------------------------------------
GH2 <- statmod::gauss.quad.prob(3, 'normal')
w <- GH2$we; v <- GH2$no

Eb2 <- mapply(function(mu, S){ # This should simply === b.hat!
  sig <- sqrt(diag(S))
  out <- matrix(NA, nr = 3, nc = 2)
  for(l in 1:3)
    out[l,] <- w[l] * (sig * v[l] + mu)
  colSums(out)
}, mu = b.hat, S = Sigma, SIMPLIFY = F)


eta <- mapply(function(X, Z, b) X[[1]] %*% beta + Z[[1]] %*% b, X = X, Z = Z, b = b.hat, SIMPLIFY = F)
tau <- mapply(function(Z, S) sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]]))), Z = Z, S = Sigma, SIMPLIFY = F)
Eexpeta2 <- mapply(function(eta, tau){
  out <- matrix(NA, nr = length(eta), nc = 3)
  for(l in 1:3)
    out[, l] <- w[l] * exp(tau * v[l] + eta)
  rowSums(out)
}, eta = eta, tau = tau, SIMPLIFY = F)


# gauss.quad.prob2 --------------------------------------------------------
# Very expensive.
Eexpeta3 <- mapply(function(eta, tau){
  out <- numeric(length(eta))
  for(q in 1:length(eta)){
    .gh <- statmod::gauss.quad.prob(3, 'normal', mu = eta[q], sigma = tau[q])
    w <- .gh$w; v <- .gh$n
    outq <- numeric(3)
    for(l in 1:3) outq[l] <- w[l] * exp(v[l])
    out[q] <- sum(outq)
  }
  out
}, eta = eta, tau = tau, SIMPLIFY = F)

# Check these three are the same ------------------------------------------
for(i in 1:250){
  exps <- cbind(Eexpeta[[i]], Eexpeta2[[i]], Eexpeta3[[i]])
  diffs <- max(abs(apply(exps, 1, diff)))
  cat(sprintf("id: %d, max difference: %.4f\n", i, diffs))
}
# Yep!


# Taking wrt joint density ------------------------------------------------
# First value only
xx <- X[[1]][[1]][1,,drop=F]; zz <- Z[[1]][[1]][1,,drop=F]; yy <- Y[[1]][[1]][1]; bb <- b.hat[[1]]

ff <- function(x, y){
  exp(xx %*% beta + zz %*% c(x,y)) * dpois(yy, exp(xx%*%beta+zz%*%b[[1]])) * 
    exp(logfti(b[[1]], S[[1]], SS[[1]][1,,drop=F], Fi[[1]],Fu[[1]][1,,drop=F],l0i=l0i[[1]],haz = l0u[[1]][1], 
           Delta = Delta[[1]], gamma_rep = gamma.rep, zeta = zeta)) * 
    mvtnorm::dmvnorm(b[[1]], rep(0,2), D)
}
pracma::integral2(Vectorize(ff), -3.118869, 3.118869, -3.221811, 3.221811)


