n <- 100

.sim <- function(ntms, theta, family){
  if(family == "gaussian"){
    D <- diag(c(0.25, 0.09))
    random.formula <- NULL
  }else if(family == "binomial"){
    D <- matrix(.40, 1, 1)
    random.formula <- list(~1)
  }else{
    D <- diag(c(0.15, 0.02))
    random.formula <- NULL
  }
  a <- simData(n = n, ntms = ntms, theta = theta,
               beta = t(c(2, -0.1, 0.1, -0.2)),
               sigma = c(0.16),
               D =  D,
               zeta = c(0, -0.2),
               family = as.list(family),
               random.formula = random.formula,
               gamma = 0.5,
               return.ranefs = TRUE)
  list(a$data, a$ranefs)
}

getSigma <- function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D){
  solve(optim(b, joint_density, joint_density_ddb,
              Y = Y, X = X, Z = Z, beta = c(2, -0.1, 0.1, -0.2), D = D, sigma = list(0.16),
              family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
              gamma_rep = rep(0.5, ncol(Z[[1]])), zeta = -0.2, beta_inds = list(0:3), b_inds = b.inds,
              K = 1L, method = 'BFGS', hessian = T)$hessian)
}

theta70ish <- c(-1.5, 0.1)
theta50ish <- c(-2, 0.1)
to.sim <- expand.grid(mi = c(5, 10, 15),
                      theta = c('medium', 'high'),
                      family = 'binomial',
                      stringsAsFactors = FALSE)

# Simulate a set of data.
sim.sets <- apply(to.sim, 1, function(x){
  m <- as.numeric(x[1]); cl <- x[2]; fam <- x[3]
  if(cl == 'low'){
    theta <- theta20ish
  }else if(cl == 'medium'){
    theta <- theta50ish
  }else{
    theta <- theta70ish
  }
  
  replicate(1, .sim(m, theta, fam), simplify = F)
})
names(sim.sets) <- sapply(apply(to.sim, 1, paste0, collapse=','), trimws)

library(Rcpp); library(RcppArmadillo)
sourceCpp('dev/metrop.cpp')

Sample.and.tau <- function(data, btrue, theta, family){
  # Longit.
  X <- Y <- Z <- vector('list', n)
  for(i in 1:n){
    X[[i]] <- Y[[i]] <- Z[[i]] <- list()
    X[[i]][[1]] <- model.matrix(~time+cont+bin, data[data$id==i,,drop=F])
    if(family!='binomial') 
      Z[[i]][[1]] <- model.matrix(~time, data[data$id==i,,drop=F])
    else
      Z[[i]][[1]] <- model.matrix(~1, data[data$id==i,,drop=F])
    Y[[i]][[1]] <- data[data$id==i,'Y.1']
  }
  b <- lapply(1:n, function(x) btrue[x,,drop=F])
  # Survival part
  fts <- sort(unique(data[data$status==1,'survtime']))
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data, center = F)
  l0 <- exp(theta[1] + theta[2] * fts)
  if(family == "binomial"){
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1|id)), parseFormula), l0)
    D <- matrix(.40, 1, 1)
    b.inds <- list(0)
    gamma.rep <- 0.5
  }else if(family == "gaussian"){
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
    D <- diag(c(0.25, 0.09))
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }else{
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
    D <- diag(c(0.15, 0.02))
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }
  
  l0u.noSurv <- lapply(sv$l0u, function(x) rep(0, length(x)))
  Del.noSurv <- as.list(rep(0, n))
  
  Omega <- list(D = D, beta = c(2, -0.1, 0.1, -0.2), sigma = list(0.16),
                gamma = 0.5, zeta = -0.2)
  
  Sigma.full <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
  },  b = b, Y = Y, X = X, Z = Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = sv$l0u)
  
  Sigma.long <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
  },  b = b, Y = Y, X = X, Z = Z, Delta = Del.noSurv, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = l0u.noSurv)
  
  tune <- if(family == 'gaussian') 6 else if(family == 'poisson') 6 else 23
  
  full <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, Sigma){
    metropolis(b, Omega, Y, X, Z, list(family), Delta, S, Fi, l0i, SS, Fu, l0u, gamma.rep,
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 500, 3500, Sigma, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = sv$l0u, Sigma = Sigma.full)
  
  FullAcc <- unlist(lapply(full, function(x) x$Acc))
  FullWalks <- lapply(full, function(x) t(x$walks))
  
  long <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, Sigma){
    metropolis(b, Omega, Y, X, Z, list(family), Delta, S, Fi, l0i, SS, Fu, l0u, gamma.rep * 0,
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 500, 3500, Sigma, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Del.noSurv, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = l0u.noSurv, Sigma = Sigma.long)
  
  LongAcc <- unlist(lapply(long, function(x) x$Acc))
  LongWalks <- lapply(long, function(x) t(x$walks))
  
  bhatLong = lapply(LongWalks, mean)
  bhatFull = lapply(FullWalks, mean)
  
  tau.full <- Map(function(Z, S) sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]]))), Z = Z, S = Sigma.full)
  tau.full.surv <- Map(function(Fu, S) sqrt(diag(tcrossprod(.5^2 * Fu %*% S, Fu))), Fu = sv$Fu, S = Sigma.full)
  tau.long <- Map(function(Z, S) sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]]))), Z = Z, S = Sigma.long)
  mu.long <- Map(function(X, Z, b) X[[1]] %*% c(2, -0.1, 0.1, -0.2) + Z[[1]] %*% b, X = X, Z = Z, b = bhatLong)
  mu.surv <- Map(function(SS, Fu, b) SS %*% Omega$zeta + Fu %*% (gamma.rep * b), SS = sv$SS, Fu = sv$Fu, b = bhatFull)
  
  list(
    FullWalks = FullWalks, FullAcc = FullAcc,
    LongWalks = LongWalks, LongAcc = LongAcc,
    bhatLong = bhatLong,
    bhatFull = bhatFull,
    tau.full = tau.full,
    tau.full.surv=tau.full.surv,
    tau.long = tau.long,
    mu.long = mu.long,
    mu.surv = mu.surv,
    sv = sv,
    surv = surv
  )
}

X <- sim.sets$`15,high,binomial`[[1]]
dat <- X[[1]]
btrue <- X[[2]]

test <- Sample.and.tau(dat, btrue, theta50ish, 'binomial')

# In EM, we only consider longitudinal part and survival densities as 
# they correspond to log-likelihoods whose expectation we need to calculate.
# Longitudinal bit only --->
GH <- gauss.quad.prob(9, 'normal')
plot(density(test$LongWalks[[1]]),  # f(b_i|Y_i;Omega{TRUE}).
     main = expression(f*"("*Y[i]*"|"*b[i]^{"TRUE"}*"; "*Omega^{"TRUE"}*")"),
     xlab = '') 
# All nodes are same over time since just random intercept -->
# points(test$tau.long[[1]][1] * GH$n, GH$w, col = 'red', pch = 'x')
# Now 'scaling up' by b.hat
for(j in 1:1){#length(test$mu.long[[1]])){
  # lines(btrue[1,] + test$tau.full[[1]][j] * GH$n, GH$w, col = 'brown', pch = 'x')
  # points(btrue[1,] + test$tau.full[[1]][j] * GH$n, GH$w, col = 'red', pch = 'x')
  # lines(btrue[1,] + test$tau.long[[1]][j] * GH$n, GH$w, col = 'magenta', pch = 'x')
  # points(btrue[1,] + test$tau.long[[1]][j] * GH$n, GH$w, col = 'blue', pch = 'x')
  lines(test$bhatFull[[1]] + test$tau.full[[1]][j] * GH$n, GH$w, col = 'brown', pch = 'x')
  points(test$bhatFull[[1]] + test$tau.full[[1]][j] * GH$n, GH$w, col = 'red', pch = 'x')
  lines(test$bhatLong[[1]] + test$tau.long[[1]][j] * GH$n, GH$w, col = 'magenta', pch = 'x')
  points(test$bhatLong[[1]] + test$tau.long[[1]][j] * GH$n, GH$w, col = 'blue', pch = 'x')
}
legend('topleft', bty = 'n', col = c( 'magenta', 'brown'), lty = c(1, 1), lwd = c(1.5,1.5),
       legend = c(
         expression("Quadrature scaled by "*hat(Sigma)[i]^{"(long)"}),
         expression("Quadrature scaled by "*hat(Sigma)[i]^{"(full)"})
       ))


# Continue tomorrow.
a <- test$FullWalks[[1]]    # f(b_i|Y_i,T_i,Delta_i;Omega{TRUE}).
plot(density(a))
lines(test$bhatFull[[1]][1] + test$tau.full.surv[[1]][1] * GH$nodes, GH$weights, col = "brown")
points(test$bhatFull[[1]][1] + test$tau.full.surv[[1]][1] * GH$nodes, GH$weights, col = "red", pch = 'x')
