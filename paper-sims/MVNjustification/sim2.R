n <- 500

.sim <- function(ntms, theta, family){
  if(family == "gaussian"){
    D <- diag(c(0.25, 0.09))
    random.formula <- NULL
  }else if(family == "binomial"){
    D <- matrix(.5, 1, 1)
    random.formula <- list(~1)
  }else{
    D <- diag(c(0.5, 0.09))
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
  
  replicate(1, .sim(m, theta, fam), simplify = F)
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
    D <- matrix(.5, 1, 1)
    b.inds <- list(0)
    gamma.rep <- 0.5
  }else if(family == "gaussian"){
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
    D <- diag(c(0.25, 0.09))
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }else{
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
    D <- diag(c(0.50, 0.09))
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
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 500, 3500, Sigma, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = sv$l0u)
  
  FullAcc <- unlist(lapply(full, function(x) x$Acc))
  FullWalks <- do.call(rbind, lapply(full, function(x) t(x$walks)))
  
  long <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    Sigma <- getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
    metropolis(b, Omega, Y, X, Z, list(family), Delta, S, Fi, l0i, SS, Fu, l0u, gamma.rep * 0,
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 500, 3500, Sigma, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Del.noSurv, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, 
  Fu = sv$Fu, l0u = l0u.noSurv)
  
  LongAcc <- unlist(lapply(long, function(x) x$Acc))
  LongWalks <- do.call(rbind, lapply(long, function(x) t(x$walks)))
  
  list(
    FullWalks = FullWalks, FullAcc = FullAcc,
    LongWalks = LongWalks, LongAcc = LongAcc
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
    ddL <- density(.SAMPLE$long[,j], na.rm = TRUE); ddLx <- ddL$x; ddLy <- ddL$y
    ddF <- density(.SAMPLE$full[,j], na.rm = TRUE); ddFx <- ddF$x; ddFy <- ddF$y
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
