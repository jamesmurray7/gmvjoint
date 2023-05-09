# #################################################################################################
# create.R                                                  Author: James Murray                  #
# Creating simulation sets (N = 100),  under a mixed trivariate joint model specification, this   #
# is achieved using simData from gmvjoint. Light-to-moderate correlation between the K={3,5}      #
# responses is induced by manually setting arguments in vech(D). Gaussian and Poisson responses   #
# are simulated under an intercept-and-slope random effects specification, where binomial ones    #
# with a random intercept only.                                                                   #
# We allow the maximal profile length to vary, along with failure rate. The number of subjects    #
# simulated is held constant at n = 250.                                                          #
# #################################################################################################

save.dir <- '~/Downloads/'              # Save to downloads (obviously change).

# Trivariate --------------------------------------------------------------
D <- diag(c(0.25, 0.09, 0.50, 0.10, 2.00))
D[1,3] <- D[3,1] <- D[1,5] <- D[5,1] <- D[3,5] <- D[5,3] <- 0.25
# check it's positive-definite?
all(eigen(D)$val > 0)
det(D) > 0
isSymmetric(D)
# Inspect correlation
round(cov2cor(D), 4)

# Parameters 
beta <- rbind(                         # Fixed effects
  c(2, -0.1, 0.1, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(1, -1, 1, -1)
)
zeta <- c(0, -0.2)                     # Time invariant survival
gamma <- c(0.5, -0.5, 0.5)             # Association
sigma <- c(0.16, 0, 0)                 # Dispersion (variance for Gaussian response)
family <- list("gaussian",
               "poisson",
               "binomial")

# Set out things that vary across simulations.
N <- 100
to.sim <- expand.grid(n = c(250), mi = c(5, 10, 15), failure = c('10%','30%', '50%'))
nms <- apply(to.sim, 1, function(x) paste0("n = ", as.numeric(x[1]), ", mi = ", as.numeric(x[2]), ", failure = ", x[3]))
sim.sets <- setNames(apply(to.sim, 1, function(x){
  n <- as.numeric(x[1]); mi <- as.numeric(x[2]); failure <- x[3]
  if(failure == '30%'){
    theta <- c(-4.,.1)    # appx. 30%
  }else if(failure == "10%"){
    theta <- c(-5.6, .1) 
  }else{
    theta <- c(-3, 0.1) # 48-54%ish
  }
  replicate(N, simData(n = n, ntms = mi, family = family, sigma = sigma, beta = beta, fup = 10,
                       D = D, gamma = gamma, zeta = zeta, theta = theta, unif.times = FALSE,
                       random.formula = list(~time, ~time, ~1))$data, simplify = F)
}), nms)

# Check average failure rate
rates <- sapply(sim.sets, function(X){
  sapply(X, function(y){
    sum(tapply(y$status, y$id, head, 1))/to.sim$n[1]
  })
})
round(colMeans(rates), 2) # On average, we get what we want!
# n = 250, mi = 5, failure = 10% n = 250, mi = 10, failure = 10% n = 250, mi = 15, failure = 10%  n = 250, mi = 5, failure = 30% 
#   0.11                            0.12                            0.11                            0.31 
# n = 250, mi = 10, failure = 30% n = 250, mi = 15, failure = 30%  n = 250, mi = 5, failure = 50% n = 250, mi = 10, failure = 50% 
#   0.31                            0.31                            0.50                            0.50 
# n = 250, mi = 15, failure = 50% 
# 0.50 

save(sim.sets, file = paste0(save.dir, 'simsets_', gsub('\\s|\\:','_',.Internal(date())), '.RData')) # About 25MB

# Five-variate ------------------------------------------------------------
rm(list=ls())
save.dir <- '~/Downloads/'     
D <- diag(c(0.25, 0.09, 0.30, 0.06,       # Gaussian (1, 2)
            0.20, 0.04, 0.50, 0.09,       # Count    (1, 2)
            2.00))                        # Binary   (1)
ints <- expand.grid(c(1,3,5,7,9), c(1,3,5,7,9))
ints <- ints[ints$Var1!=ints$Var2,]
D[cbind(ints$Var1, ints$Var2)] <- 0.125
all(eigen(D)$val > 0) && det(D) > 0 && isSymmetric(D)
# Inspect correlation
round(cov2cor(D), 4)

# Parameters 
beta <- rbind(                            # Fixed effects ----
  c(2, -0.1, 0.1, -0.2),                  # Gaussian (1)
  -c(2, -0.1, 0.1, -0.2),                 # Gaussian (2)
  c(2, -0.1, 0.1, -0.2),                  # Count (1)
  c(2, -0.1, 0.1, -0.2),                  # Count (2)
  c(1, -1, 1, -1)                         # Binomial (1)
)

zeta <- c(0, -0.2)                        # Time invariant survival
gamma <- c(0.25, -0.25, 0.25, -0.25, 0.30)# Association
sigma <- c(0.16, 0.16, 0, 0, 0)           # Dispersion (variance for Gaussian response)
family <- list("gaussian", "gaussian",
               "poisson", "poisson",
               "binomial")
random.formula <- list(~time, ~time, ~time, ~time, ~1)

sim.sets <- replicate(500,
                      simData(family = family, sigma = sigma, beta = beta, D = D, gamma = gamma, 
                              zeta = zeta, theta = c(-2.900088, 0.1),
                              random.formula = random.formula, n = 500, ntms = 10)$data,
                      simplify = F)

save(sim.sets, file = paste0(save.dir, 'simsets_5variate_', gsub('\\s|\\:','_',.Internal(date())), '.RData')) # About 50MB

# fit <- joint(list(
#   Y.1 ~ time + cont + bin + (1 + time|id),
#   Y.2 ~ time + cont + bin + (1 + time|id),
#   Y.3 ~ time + cont + bin + (1 + time|id),
#   Y.4 ~ time + cont + bin + (1 + time|id),
#   Y.5 ~ time + cont + bin + (1|id)
# ), Surv(survtime, status) ~ bin, data, family  = family, 
# control = list(verbose=T))
