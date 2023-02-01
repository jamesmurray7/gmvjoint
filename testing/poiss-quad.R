rm(list=ls())
# Simulate some data ------------------------------------------------------

n <- 100; mi <- 5
time <- rep(0:(mi-1), n)
D <- diag(c(.25, .04))
b <- MASS::mvrnorm(100, c(0,0), D)
df <- data.frame(id = rep(1:n, each = mi),
                 time = time,
                 cont = rep(rnorm(n), each = mi),
                 bin = rep(rbinom(n,1,.5),each = mi))
head(df)
X <- model.matrix(~time+cont+bin, df)
Z <- model.matrix(~time,df)

eta <- X %*% c(1,0.2,0.4,-0.2) + rowSums(Z * b[df$id,])
Y <- rpois(n*mi, exp(eta))

df$Y <- Y
head(df)

# f(y, b)
true.ll <- sum(dpois(Y, exp(eta), T)) + sum(mvtnorm::dmvnorm(b, rep(0,2), D, log = T)) 

library(glmmTMB)
mod <- glmmTMB(Y~time+cont+bin+(1+time|id),df,poisson)
mod.ll <- c(logLik(mod))
beta <- fixef(mod)$cond
bmod <- as.matrix(ranef(mod)$cond$id)
Dmod <- VarCorr(mod)$cond$id
etamod <- X %*% beta + rowSums(Z * bmod[df$id,])

tosplit <- split(1:(n*mi),df$id)
Ylist <- lapply(tosplit, function(x) Y[x])
Zlist <- lapply(tosplit, function(x) Z[x,])
Xlist <- lapply(tosplit, function(x) X[x,])
blist <- lapply(1:n, function(x) bmod[x,,drop=F])
etamodlist <- lapply(tosplit, function(x) etamod[x,,drop=F])

# GH
ngh <- 3
GHfun <- function(ngh, D){
  # Xform D into sigmas (SDs) of normal distn
  sigs <- sqrt(diag(D)); q <- length(sigs)
  nodes <- vector('list', q)
  for(l in 1:q){
    GH <- statmod::gauss.quad.prob(ngh, 'normal', sigma = sigs[l])
    nodes[[l]] <- GH$nodes
    weights <- GH$weights
  }
  list(nodes = do.call(cbind, nodes), weights = weights)
}

# function to define f(y,b) (= f(y|b) * f(b))
f.yb <- function(b, X, Y, Z, beta, D){
  eta <- X %*% beta + Z %*% c(b)
  pois <- dpois(Y, exp(eta), FALSE)
  dnorm <- mvtnorm::dmvnorm(b, c(0,0), D, F)
  prod(pois * dnorm)
}
(mapply(function(b, X, Y, Z) f.yb(b, X, Y, Z, beta, D), 
       b = blist, X = Xlist, Y = Ylist, Z = Zlist))

GH.appx <- mapply(function(b, X, Y, Z){
  gauher <- GHfun(ngh, Dmod)
  w <- gauher$w; v <- gauher$no
  fyb <- numeric(ngh)
  for(l in 1:ngh) fyb[l] <- w[l] * 
    f.yb(b = v[l,,drop=F], X = X, Y = Y, Z = Z, beta = beta, D = Dmod)
  sum(fyb)
}, b = blist, X = Xlist, Y = Ylist, Z = Zlist)

sum(log(GH.appx)) # Reasonable appx. 


# Adaptive GH -------------------------------------------------------------
# Find modes and curvatures of random effects
logf.yb <- function(b, X, Y, Z, beta, D){
  eta <- X %*% beta + Z %*% c(b)
  -sum(dpois(Y, exp(eta), T) + mvtnorm::dmvnorm(b, c(0,0), D, T))
}

modes.H <- mapply(function(b, X, Y, Z){
  opt <- optim(b, logf.yb, X = X, Y = Y, Z = Z, beta = beta,
               D = Dmod, method = 'BFGS', hessian = TRUE)
  list(mode = opt$par, Hess = opt$hess)
}, b = blist, X = Xlist, Y = Ylist, Z = Zlist, SIMPLIFY = F)

# mode
bhat <- lapply(modes.H, el)
# Lower triangle Cholesky
Ch <- lapply(modes.H, function(x) chol(x$Hess))

AGH.appx <- mapply(function(b, X, Y, Z, C){
  GH <- GHfun(ngh, solve(C)^2)
  w <- GH$w; v <- GH$n
  fyb <- numeric(ngh)
  for(l in 1:ngh) fyb[l] <- w[l] * f.yb(b = b + v[l,,drop=T], X = X, Y = Y, Z = Z, beta = beta, D = Dmod)
  # for(l in 1:ngh) fyb[l] <- w[l] * f.yb(b = v[l,,drop=F] + b, X = X, Y = Y, Z = Z, beta = beta, D = Dmod)
  exp(determinant.matrix(C,T)$mod) * sum(fyb)
}, b = bhat, X = Xlist, Y = Ylist, Z = Zlist, C = Ch)

sum(log(AGH.appx))
