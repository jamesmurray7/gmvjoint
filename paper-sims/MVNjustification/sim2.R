n <- 250

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

# Some good thetas to use.
theta70ish <- c(-1.5, 0.1)
theta50ish <- c(-2, 0.1)
theta20ish <- c(-3, 0.1)
to.sim <- expand.grid(mi = c(5, 10, 15),
                      theta = c('low', 'medium', 'high'),
                      family = c('gaussian', 'poisson', 'binomial'),
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
  
  .sim(m, theta, fam)
})
names(sim.sets) <- sapply(apply(to.sim, 1, paste0, collapse=','), trimws)

library(Rcpp); library(RcppArmadillo)
sourceCpp('dev/metrop.cpp')

Sample <- function(data, btrue, theta, family){
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
  
  tune <- if(family == 'gaussian') 6 else if(family == 'poisson') 6 else 23
  
  full <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    Sigma <- getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
    # print(l0i)
    metropolis(b, Omega, Y, X, Z, list(family), Delta, S, Fi, l0i, SS, Fu, l0u, gamma.rep,
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 250, 2000, Sigma, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = sv$l0u)
  
  FullAcc <- unlist(lapply(full, function(x) x$Acc))
  FullWalks <- do.call(rbind, lapply(full, function(x) t(x$walks)))
  
  list(
    FullWalks = FullWalks, FullAcc = FullAcc
  )
}

# Sample from f(b|...) 
to.apply <- function(x){ # x is a ROW from to.sim
  mi <- as.numeric(x[1]); fr <- x[2]; family <- as.character(x[3])
  if(fr == 'low'){
    theta <- theta20ish
  }else if(fr == 'medium'){
    theta <- theta50ish
  }else{
    theta <- theta70ish
  }
  
  if(family == "poisson")
    D <- diag(c(0.15, 0.02))
  else if(family == "gaussian")
    D <- diag(c(0.25, 0.09))
  else
    D <- matrix(.40,1,1)
  
  .sim.sets.lookup <- trimws(paste0(mi,',',fr,',',family))
  X <- sim.sets[[.sim.sets.lookup]]
  data <- X[[1]]
  btrue <- X[[2]]
  
  .SAMPLE <- Sample(data, btrue, theta, family)
  message(.sim.sets.lookup, " done!")
  
  .SAMPLE
}

# Create some plots.

library(dplyr)
library(tidyr)
library(ggplot2)
source('zzz/theme_csda.R')

# Consider each family in turn, and DELETE at end of each, 
# otherwise I run out of RAM!

# Function to create a panelled plot
plotter <- function(f){
  file.name <- paste0(f,'.png')
  inds <- which(to.sim$family == f)
  this.to.plot <- to.sim[inds,]
  samples <- apply(to.sim[inds,], 1, to.apply)
  
  samples <- do.call(rbind, lapply(seq_along(samples), function(i){
    W <- samples[[i]]$FullWalks
    if(f == 'binomial')
      colnames(W) <- 'b[0]'
    else
      colnames(W) <- c('b[0]', 'b[1]')
    W <- as.data.frame(W)
    W$r <- this.to.plot$mi[i]
    W$o <- this.to.plot$theta[i]
    W
  }))
  
  gg.object <- samples %>% 
    mutate(r = factor(r, c('5','10','15')),
           omega.fac = case_when(
             o == 'low' ~ 'omega == 0.2',
             o == 'medium' ~ 'omega == 0.5',
             o == "high" ~ 'omega == 0.7'
    )) %>% 
    mutate(omega.fac = factor(omega.fac, c('omega == 0.2', 
                                           'omega == 0.5',
                                           'omega == 0.7'))) %>% as_tibble() %>% 
    select(-o) %>% 
    pivot_longer(-c(`r`, `omega.fac`)) %>% 
    ggplot(aes(x = value, colour = r)) + 
    stat_density(geom = 'line', position = 'identity') + 
    scale_colour_brewer(palette = 'Dark2') + 
    theme_csda() + 
    # facet_grid(omega.fac~name, labeller = label_parsed) + 
    facet_wrap(omega.fac~name, ncol = 2, labeller = label_parsed, scales = 'free') +
    labs(y = 'Density', x = '', colour = expression(r~'=')) 
    
  gg.object
  # save(gg.object, file = paste0('./paper-sims/MVNjustification/',file.name,'.RData'))
  ggsave(paste0('./paper-sims/MVNjustification/',file.name), width = 140, height = 90, units = 'mm')
    
  # Remove to free up RAM.
  rm(samples)
  dev.off()
  on.exit(gc())
}

plotter('poisson')
plotter('binomial')
