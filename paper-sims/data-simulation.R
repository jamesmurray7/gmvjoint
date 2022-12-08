rm(list=ls())
saveDir <- paste0(getwd(),'/Simulations/data/')
source('simData.R')


# Simulation of K={1,2,3} Gaussian responses ------------------------------
K <- c(1, 2, 3)
beta <- do.call(rbind, replicate(3, c(0.2, -0.1, 0.1, -0.2), simplify = F))
gamma <- c(0.5, -0.5, 0.5)
var.e <- 0.16

D <- diag(c(0.16, 0.16, 0.25, 0.25, 0.25, 0.16))
ltriD <- c(0.08, 0.02, 0.04, 0.00, -0.04, 
           0.04, 0.00, -0.08, -0.08, 
           0.10, 0.05, 0.02, 
           0.05, -0.04,
           0.10)
D[lower.tri(D, F)] <- ltriD
D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
isSymmetric(D); any(eigen(D)$val<0);det(D)<=0

zeta <- c(0, -0.20) # drug effect only!
inds <- split(1:6, rep(c(1,2,3), each = 2))
for(k in K){
  # beta and family
  beta.k <- beta[1:k,,drop=F]
  family.k <- as.list(rep('gaussian', k))
  # pick-out submatrix of D
  inds.k <- do.call(c, lapply(1:k, function(x) inds[[x]]))
  D.k <- D[inds.k, inds.k]
  cat(paste('k =', k))
  simulated.data <- replicate(100, 
                              simData(n = 250, ntms = 15, family = family.k, 
                                      beta = beta.k, var.e = var.e, D = D.k,
                                      gamma = gamma[1:k], zeta = zeta,
                                      theta = c(-3.25, 0.2)
                                      ),
                              simplify = F)
  cat('\n\n')
  save(simulated.data, file = paste0(saveDir, 'gaussianK-', k, '.RData'))
}


# Simulation of K={1,2,3} Count responses ---------------------------------
K <- 1:3
beta <- do.call(rbind, replicate(3, c(2, -0.1, 0.1, -0.2), simplify = F))
gamma <- c(0.3, -0.3, 0.3)
D <- diag(c(0.16, 0.09, 0.25, 0.16, 0.25, 0.16))
ltriD <- c(0.06, 0.02, 0.04, 0.00, -0.04, 
           0.03, 0.00, -0.06, 0.00,
           0.08, 0.05, 0.04, 
           0.04, -0.04, 
           0.12)
D[lower.tri(D, F)] <- ltriD
D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
isSymmetric(D); any(eigen(D)$val<0);det(D)<=0
zeta <- c(0.0, -0.2)
inds <- split(1:6, rep(c(1,2,3), each = 2))
for(k in K){
  # beta and family
  beta.k <- beta[1:k,,drop=F]
  family.k <- as.list(rep('poisson', k))
  # pick-out submatrix of D
  inds.k <- do.call(c, lapply(1:k, function(x) inds[[x]]))
  D.k <- D[inds.k, inds.k]
  cat(paste('k =', k, '\n'))
  simulated.data <- replicate(100, 
                              simData(n = 250, ntms = 15, family = family.k, 
                                      beta = beta.k, var.e = var.e, D = D.k,
                                      gamma = gamma[1:k], zeta = zeta,
                                      theta = c(-3.10, 0.2)
                              ),
                              simplify = F)
  cat('\n')
  save(simulated.data, file = paste0(saveDir, 'poissonK-', k, '.RData'))
}

# Simulation of K={1,2,3} Count responses ---------------------------------
K <- 1:3 ### MIGHT HAVE TO DO THIS WITH NO RANDOM SLOPE!
beta <- do.call(rbind, replicate(3, c(1, -0.1, 0.1, -0.2), simplify = F))
gamma <- c(0.3, -0.3, 0.3)
D <- diag(c(0.16, 0.09, 0.25, 0.16, 0.25, 0.16))
ltriD <- c(0.06, 0.02, 0.04, 0.00, -0.04, 
           0.03, 0.00, -0.06, 0.00,
           0.08, 0.05, 0.04, 
           0.04, -0.04, 
           0.12)
D[lower.tri(D, F)] <- ltriD
D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
isSymmetric(D); any(eigen(D)$val<0);det(D)<=0
zeta <- c(0.0, -0.2)
inds <- split(1:6, rep(c(1,2,3), each = 2))
for(k in K){
  # beta and family
  beta.k <- beta[1:k,,drop=F]
  family.k <- as.list(rep('binomial', k))
  # pick-out submatrix of D
  inds.k <- do.call(c, lapply(1:k, function(x) inds[[x]]))
  D.k <- D[inds.k, inds.k]
  cat(paste('k =', k, '\n'))
  simulated.data <- replicate(100, 
                              simData(n = 250, ntms = 15, family = family.k, 
                                      beta = beta.k, var.e = var.e, D = D.k,
                                      gamma = gamma[1:k], zeta = zeta,
                                      theta = c(-3.10, 0.2)
                              ),
                              simplify = F)
  cat('\n')
  save(simulated.data, file = paste0(saveDir, 'binomialK-', k, '.RData'))
}
