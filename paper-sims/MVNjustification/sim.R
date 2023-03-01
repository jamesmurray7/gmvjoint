.sim <- function(ntms, theta, family){
  a <- simData(n = 250, ntms = ntms, theta = theta,
               beta = t(c(2, -0.1, 0.1, -0.2)),
               D =  diag(c(0.25, 0.05)),
               zeta = c(0, -0.2),
               family = as.list(family),
               gamma = 0.5,
               return.ranefs = TRUE)
  list(a$data, a$ranefs)
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
  
  replicate(1, .sim(m, theta, fam), simplify = F)
})
names(sim.sets) <- sapply(apply(to.sim, 1, paste0, collapse=','), trimws)

library(MCMCpack)
dens <- function(b, Y, X, Z, S, SS, Fu, Fi, l0i, haz, Del, gam, zet,
                 family,survival = TRUE){
  ff <- family[[1]]
  eta <- X[[1]] %*% c(2, -0.1, 0.1, -0.2) + Z[[1]] %*% b
  longitdens <- switch(ff,
                       gaussian = prod(dnorm(Y[[1]], eta, sqrt(0.16))),
                       poisson = prod(dpois(Y[[1]], exp(eta))),
                       binomial = prod(dbinom(Y[[1]], 1, plogis(eta))))
  REdens <- mvtnorm::dmvnorm(b, mean = rep(0, 2), sigma = diag(c(0.25, 0.05)))
  out <- longitdens * REdens
  if(survival) out <- out * exp(logfti(b, S, SS, Fi, Fu, l0i, haz, Del, gam, zet))
  out
}

Sample <- function(data, theta, family){
  # Longit.
  X <- Y <- Z <- vector('list', 250)
  for(i in 1:250){
    X[[i]] <- Y[[i]] <- Z[[i]] <- list()
    X[[i]][[1]] <- model.matrix(~time+cont+bin, data[data$id==i,,drop=F])
    Z[[i]][[1]] <- model.matrix(~time, data[data$id==i,,drop=F])
    Y[[i]][[1]] <- data[data$id==i,'Y.1']
  }
  
  # Survival part
  fts <- sort(unique(data[data$status==1,'survtime']))
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data, center = F)
  l0 <- exp(theta[1] + theta[2] * fts)
  sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)
  
  full <- Map(function(X, Y, Z, S, SS, Fu, Fi, l0i, haz, Del){
    a <- tryCatch(MCMCmetrop1R(dens, theta.init = c(0,0), 
                               Y = Y, X = X, Z = Z,
                               S = S, SS = SS, Fu = Fu, Fi = Fi, l0i = l0i,
                               haz = haz, Del = Del, gam = rep(0.5, 2), zet = -0.2, survival = TRUE,
                               family = family,
                               burnin = 50, mcmc = 3500, logfun = F, tune = 2.5),
                  error = function(e) NA
    )
    a
  }, X = X, Y = Y, Z = Z, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0i = sv$l0i,
  haz = sv$l0u, Del = surv$Delta)
  
  long <- Map(function(X, Y, Z, S, SS, Fu, Fi, l0i, haz, Del){
    a <- tryCatch(MCMCmetrop1R(dens, theta.init = c(0,0), 
                               Y = Y, X = X, Z = Z,
                               S = S, SS = SS, Fu = Fu, Fi = Fi, l0i = l0i,
                               haz = haz, Del = Del, gam = rep(0.5, 2), zet = -0.2, survival = FALSE,
                               family = family,
                               burnin = 50, mcmc = 3500, logfun = F, tune = 2.5),
                  error = function(e) NA
    )
    a
  }, X = X, Y = Y, Z = Z, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0i = sv$l0i,
  haz = sv$l0u, Del = surv$Delta)
  
  list(long = do.call(rbind, long), full = do.call(rbind, full))
  
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
  
  .sim.sets.lookup <- trimws(paste0(mi,',',fr,',',family))
  X <- sim.sets[[.sim.sets.lookup]][[1]]
  data <- X[[1]]
  
  .SAMPLE <- Sample(data, theta, family)
  message(.sim.sets.lookup, " done!")
  
  # Save an image of the plot
  png(filename = paste0('paper-sims/MVNjustification/output/',
                        gsub('\\,','_',.sim.sets.lookup), '.png'),
      width = 140, height = 75, units = 'mm', pointsize = 9, res = 2000)
  par(mfrow = c(1,2))
  for(j in 1:2){
    ddL <- density(.SAMPLE$long[,j]); ddLx <- ddL$x; ddLy <- ddL$y
    ddF <- density(.SAMPLE$full[,j]); ddFx <- ddF$x; ddFy <- ddF$y
    xs <- seq(min(pmin(ddLx, ddFx)), max(pmax(ddLx, ddFx)), length.out = 1000)
    DN <- dnorm(xs, 0, sqrt(c(.25, .05)[j]))
    ylims <- c(min(pmin(ddLy, ddFy), min(DN)), max(pmax(ddLy, ddFy), max(DN)))
    plot(ddL, xlim = c(min(xs), max(xs)), ylim = ylims,
         xlab = bquote(b[.(j-1)]), main = '')
    lines(ddF, col = 'steelblue')
    lines(DN~xs, lty = 3, col = 'red3')
  }
  dev.off()
  
  .SAMPLE
}


# Apply over ROWS of to.sim...
# Run overnight..!
ALL <- apply(to.sim, 1, to.apply)
