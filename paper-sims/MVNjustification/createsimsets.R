all.sims <- expand.grid(mi = c(6,10,15), 
                        n =  c(100, 250, 500))
# Families
family <- list("gaussian", "poisson", "binomial")

# VarCorr matrix for REs
D <- diag(c(0.16, 0.04, 0.25, 0.05, 0.25))
D[lower.tri(D, )] <- c(0.03, 0.02, 0.04, 0.0,  # Gaussian (Intercept) and others.
                       0.03, 0.00, -0.06,      # Gaussian (slope) and others.
                       0.08, 0.05,             # Poisson (Intercept) and others.
                       0.04)                   # Poisson slope and binomial intercept

D[upper.tri(D)] <- t(D)[upper.tri(D)]

# Check this ok
isSymmetric(D) && all(eigen(D)$val>0) && det(D) > 0

D <- diag(c(0.16, 0.04, 0.25, 0.05, 0.20, 0.04))
D[lower.tri(D, )] <- c(0.03, 0.02, 0.04, 0.0, 0.0, # Gaussian (Intercept) and others.
                       0.03, 0.00, -0.06,  0.0  ,  # Gaussian (slope) and others.
                       0.08, 0.05,   0.0,          # Poisson (Intercept) and others.
                       0.04, 0.0,0.0)                   # Poisson slope and binomial intercept

# Fixed effects
beta <- rbind(
  c(0.2, -0.2, 0.2, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(2, -0.1, 0.1, -0.2)
  # c(1, -1, 1, -1)
)

# Dispersion (only needed for sigma)
sigma <- c(0.16, 0, 0)
sigma <- c(0.16, 0, 0.25)

# gamma & zeta
gamma <- c(-0.5, 0.25, 0.40)
zeta <- c(-0.0, 0.3)

# Specify random effects formulae (W_k)
random.formula <- list(~time, ~time, ~1)

all.sims2 <- setNames(lapply(1:nrow(all.sims), function(i){
  n <- all.sims[i,'n']; mi <- all.sims[i, 'mi']
  replicate(100, simData(family = list('gaussian', 'poisson', 'poisson'), 
                         sigma = sigma, ntms = mi, n = n,
                         beta = beta, D = D, gamma = gamma, zeta = zeta, 
                         # random.formula = random.formula,
                         theta = c(-2.5, 0.1))$data, simplify = F)
}), apply(all.sims, 1, 
          function(x) paste0('n = ', x[2], ', mi = ', x[1])))

ff <- function(d){ # check this!
  joint(list(Y.1~time+cont+bin + (1 + time|id),
             Y.2~time+cont+bin + (1 + time|id),
             Y.3~time+cont+bin + (1+time|id)),
        Surv(survtime, status) ~ bin,
        data = d, family = list('gaussian', 'poisson', 'poisson'))
}

fits <- setNames(vector('list', nrow(all.sims)), names(all.sims2))
for(s in 1:length(fits)){
  message("\n=======\nn=",all.sims[s,'n'],' mi=',all.sims[s,'mi'],'\n=======\n')
  sset <- all.sims2[[s]] # 100-list of data
  fit <- mclapply(sset, ff)
  cat(sprintf("%d out of %d done.\r", s, length(fits)))
}

# all.sims2 <- setNames(lapply(1:nrow(all.sims), function(i){
#   n <- all.sims[i,'n']; mi <- all.sims[i, 'mi']
#   replicate(100, simData(family = family, sigma = sigma, ntms = mi, n = n,
#                          beta = beta, D = D, gamma = gamma, zeta = zeta, random.formula = random.formula,
#                          theta = c(-2.5, 0.1))$data, simplify = F)
# }), apply(all.sims, 1, 
#           function(x) paste0('n = ', x[2], ', mi = ', x[1])))

test <- setNames(lapply(1:nrow(all.sims), function(i){
  n <- all.sims[i,'n']; mi <- all.sims[i, 'mi']
  print(mi)
  if(mi == 6) theta <- c(-2.3, 0.1)
  if(mi == 10) theta <- c(-3.0, 0.1)
  if(mi == 15) theta <- c(-3.8, 0.1)
  replicate(2, simData(sigma = sigma, 
                       ntms = mi, n = n, fup = mi -1,
                       beta = beta, D = D, gamma = gamma, zeta = zeta,
                       family = list('gaussian', 'poisson', 'poisson'),
                       # family = family,
                       # random.formula = random.formula,
                      theta = theta)$data, simplify = F)
}), apply(all.sims, 1, 
          function(x) paste0('n = ', x[2], ', mi = ', x[1])))

library(parallel)
a <- lapply(test, function(x){
  out <- mclapply(x, function(d){
    joint(list(Y.1~time+cont+bin + (1 + time|id),
               Y.2~time+cont+bin + (1 + time|id),
               Y.3~time+cont+bin + (1|id)),
          Surv(survtime, status) ~ bin,
          data = d, family = list('gaussian', 'poisson', 'poisson'))
  out
}, mc.cores = 2L)})

ff <- function(d){
  joint(list(Y.1~time+cont+bin + (1 + time|id),
             Y.2~time+cont+bin + (1 + time|id),
             Y.3~time+cont+bin + (1+time|id)),
             Surv(survtime, status) ~ bin,
             data = d, family = list('gaussian', 'poisson', 'poisson'))
}

