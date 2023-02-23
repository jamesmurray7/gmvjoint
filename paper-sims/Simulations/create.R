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
vD <- c(0.25, 0.04, -0.02, 0.00, 0.03, # Gaussian intercept and G-slope; P-int+ P-slope; B-int.
        0.09, 0.00, -0.01, 0.00,       # Gaussian slope and with P-int+ P-slope; B-int.
        0.25, -0.04, 0.02,             # Poisson intercept with P-slope; B-int.
        0.05, 0.00,
        0.10)
D <- vech2mat(vD, 5)        
# check it's positive-definite?
all(eigen(D)$val > 0)
det(D) > 0
# Inspect correlation
round(cov2cor(D), 4)

# Parameters 
beta <- rbind(                         # Fixed effects
  c(2, -0.5, 0.1, 0.5),
  -c(2, -0.5, 0.1, 0.5),
  c(1, -1, 1, -1)
)
zeta <- c(0, -0.2)                     # Time invariant survival
gamma <- c(-0.5, 0.5, -0.5)            # Association
sigma <- c(0.16, 0, 0)                 # Dispersion
family <- list("gaussian",
               "poisson",
               "binomial")

# Set out things that vary across simulations.
n <- 500; N <- 100
to.sim <- expand.grid(mi = c(5, 10, 15), failure = c('low', 'medium', 'high'))
nms <- apply(to.sim, 1, function(x) paste0("mi = ", as.numeric(x[1]), ", failure = ", x[2]))
sim.sets <- setNames(apply(to.sim, 1, function(x){
  mi <- as.numeric(x[1]); failure <- x[2]
  if(failure == 'low'){
    theta <- c(-3.3, 0.1)
  }else if(failure == 'medium'){
    theta <- c(-2.2, 0.1)
  }else if(failure == 'high'){
    theta <- c(-1.5, 0.1)
  }
  replicate(N, simData(n = n, ntms = mi, family = family, sigma = sigma, beta = beta,
                       D = D, gamma = gamma, zeta = zeta, theta = theta,
                       random.formula = list(~time, ~time, ~1))$data, simplify = F)
}), nms)

save(sim.sets, file = paste0(save.dir, 'simsets_', gsub('\\s','_',.Internal(date())), '.RData')) # About 50MB

# Five-variate ------------------------------------------------------------

# Ongoing =-=-=-=-=-=-=-=-=-=

