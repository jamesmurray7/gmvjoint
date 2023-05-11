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
  uu <- optim(b, joint_density, joint_density_ddb,
        Y = Y, X = X, Z = Z, beta = c(2, -0.1, 0.1, -0.2), D = D, sigma = list(0.16),
        family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
        gamma_rep = rep(0.5, ncol(Z[[1]])), zeta = -0.2, beta_inds = list(0:3), b_inds = b.inds,
        K = 1L, method = 'BFGS', hessian = T)
  list(bhat = uu$par, Sigma = solve(uu$hessian))
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

Sample <- function(data, btrue, theta, family, ids){
  # Longit.
  X <- Y <- Z <- setNames(vector('list', length(ids)), paste0("id: ", ids))
  for(i in seq_along(ids)){
    X[[i]] <- Y[[i]] <- Z[[i]] <- list()
    X[[i]][[1]] <- model.matrix(~time+cont+bin, data[data$id==ids[i],,drop=F])
    if(family!='binomial') 
      Z[[i]][[1]] <- model.matrix(~time, data[data$id==ids[i],,drop=F])
    else
      Z[[i]][[1]] <- model.matrix(~1, data[data$id==ids[i],,drop=F])
    Y[[i]][[1]] <- data[data$id==ids[i],'Y.1']
  }
  b <- lapply(seq_along(ids), function(x) btrue[ids[x],,drop=F])
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
  
  Omega <- list(D = D, beta = c(2, -0.1, 0.1, -0.2), sigma = list(0.16),
                gamma = 0.5, zeta = -0.2)
  
  tune <- if(family == 'gaussian') 6 else if(family == 'poisson') 6 else 23
  
  Delta <- lapply(seq_along(ids), function(x) surv$Delta[[ids[x]]])
  S <- lapply(seq_along(ids), function(x) sv$S[[ids[x]]])
  Fi <- lapply(seq_along(ids), function(x) sv$Fi[[ids[x]]])
  l0i <- lapply(seq_along(ids), function(x) sv$l0i[[ids[x]]])
  SS <- lapply(seq_along(ids), function(x) sv$SS[[ids[x]]])
  Fu <- lapply(seq_along(ids), function(x) sv$Fu[[ids[x]]])
  l0u <- lapply(seq_along(ids), function(x) sv$l0u[[ids[x]]])
  
  Sigma <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    Sigma <- getSigma(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, 
  Fu = Fu, l0u = l0u)
  
  b.hat <- lapply(Sigma, el, 1)
  Sigma <- lapply(Sigma, el, 2)
  
  full <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, Sig){
    metropolis(b, Omega, Y, X, Z, list(family), Delta, S, Fi, l0i, SS, Fu, l0u, gamma.rep,
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 250, 10000, Sig, tune)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, 
  Fu = Fu, l0u = l0u, Sig = Sigma)
  
  Acc <- unlist(lapply(full, function(x) x$Acc))
  Walks <- setNames(lapply(full, function(x) t(x$walks)), names(X))
  
  list(
    Walks = Walks, b.hat = b.hat, Sigma = Sigma, Acc = Acc,
    family = family, theta = theta,
    nums = sapply(Z, function(x) nrow(x[[1]]))
  )
}

plotSample <- function(S){
  ids <- as.numeric(gsub("id\\:\\s", "", names(S$Walks)))
  
  if(S$family == "poisson")
    D <- diag(c(0.15, 0.02))
  else if(S$family == "gaussian")
    D <- diag(c(0.25, 0.09))
  else
    D <- matrix(.40,1,1)
  
  par(mfrow=c(length(ids), ncol(D)))
  for(i in seq_along(ids)){
    ii <- ids[i]
    
    for(j in 1:ncol(D)){
      dens <- density(S$Walks[[i]][,j])
      plot(dens, main = "", ylab = names(S$Walks)[i], xlab ="")
      xmin <- min(dens$x); xmax <- max(dens$x)
      curve(dnorm(x, mean = S$b.hat[[i]][j], sd = sqrt(S$Sigma[[i]][j,j])),
            from = xmin, to = xmax,
            col = 'tomato', lty = 5, add = TRUE, )
    }
  }
  
  on.exit(gc())
}


