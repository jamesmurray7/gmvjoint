theta70ish <- c(-1.5, 0.1)
theta50ish <- c(-2, 0.1)
theta20ish <- c(-3, 0.1)
for(f in dir("./paper-sims/MVNjustification/output", full.names = T)) unlink(f)
theta <- theta50ish
# Function to obtain a simulated set of univariate joint data.
.sim <- function(family){
  if(family == "gaussian"){
    D <- diag(c(0.25, 0.09))
    random.formula <- NULL
  }else if(family == "binomial"){
    D <- matrix(2, 1, 1)
    random.formula <- list(~1)
  }else{
    D <- diag(c(0.25, 0.09))
    random.formula <- NULL
  }
  a <- simData(n = 250, ntms = 15, theta = theta,
               beta = t(c(2, -0.1, 0.1, -0.2)),
               sigma = c(0.16),
               D =  D,
               zeta = c(0, -0.2),
               family = as.list(family),
               random.formula = random.formula,
               gamma = 0.5,
               return.ranefs = TRUE)
  list(a$data, a$ranefs, D)
}

# Obtain b.hat and Sigma.hat given observed data at TRUE parameter estimates.
getSigma <- function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u, family, b.inds, D){
  uu <- optim(b, joint_density, joint_density_ddb,
              Y = Y, X = X, Z = Z, beta = c(2, -0.1, 0.1, -0.2), D = D, sigma = list(0.16),
              family = as.list(family), Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u,
              gamma_rep = rep(0.5, ncol(Z[[1]])), zeta = -0.2, beta_inds = list(0:3), b_inds = b.inds,
              K = 1L, method = 'BFGS', hessian = T)
  list(bhat = uu$par, Sigma = solve(uu$hessian))
}

# Function to sample given data, true random effects, known family and target ids.
Sample <- function(data, btrue, family, ids){
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
  
  # Tuning parameters, given D, these return about 22.5% acceptane.
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
               list(0:3), b.inds, 1L, length(b.inds[[1]]), 1000, 10000, Sig, tune)
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

# Some good thetas to use.
plotFamily <- function(family, save.dir = './paper-sims/MVNjustification/output/'){
  d <- .sim(family)
  data <- d[[1]]; b.true <- d[[2]]; D <- d[[3]]
  q <- ncol(D)
  # Find the _first_ subject who has five, ten or fifteen follow-up visits.
  r <- with(data, tapply(time, id, length))
  r5 <- unname(which(r==5)[1]); r10 <- unname(which(r==10)[1]); r15 <- unname(which(r==15)[1])
  ids <- c(r5, r10, r15)
  ss <- Sample(data = data, btrue = b.true, family = family, ids = ids)
  # Create empty lists to store densities, etc.
  densinfos <- setNames(vector('list', 3), paste0("r",c(5,10,15)))
  
  # First: Create a 3 x q (for Gaussian/Poisson) or 1x3 panel for Binomial
  fn <- paste0(save.dir,'panelled_', family)
  fn <- paste0(fn, '.png')
  if(family != 'binomial')
    png(filename = fn, width = 140, height = 120, units = 'mm', pointsize = 9, res = 2000)
  else
    png(filename = fn, width = 140, height = 60, units = 'mm', pointsize = 9, res = 2000)
  if(family!="binomial") par(mfrow=c(3, q)) else par(mfrow=c(1, 3))
  for(i in seq_along(ids)){
    ii <- ids[i]
    b.hat <- ss$b.hat[[i]]
    Sigma <- ss$Sigma[[i]]
    ii.r <- unname(ss$nums[i])
    densinfos[[i]] <- list(cond.dens = setNames(vector("list", q), paste0("b",0:(q-1))),
                           dnorm.dens = setNames(vector("list", q), paste0("b",0:(q-1))),
                           ylims = setNames(vector("list", q), paste0("b",0:(q-1))),
                           xs = setNames(vector("list", q), paste0("b",0:(q-1))),
                           b.hat = b.hat,
                           Sigma = Sigma)
    
    for(j in 1:q){
      Sigma.jj <- sqrt(Sigma[j,j])
      b.hat.jj <- b.hat[j]
      dens <- density(ss$Walks[[i]][,j])
      xmin <- min(dens$x); xmax <- max(dens$x)
      xs <- seq(xmin, xmax, length.out=1e3)
      ymin <- min(dens$y); ymax <- max(dens$y)
      dn <- dnorm(dens$x, mean = b.hat.jj, sd = Sigma.jj)
      ylims <- c(0, max(ymax, max(dn)))
      plot(dens, xlab = '', ylab = bquote(m[i]==.(ii.r)), main = bquote(b[.(j-1)]),
           ylim = ylims)
      lines(dn ~ dens$x, lty = 5, col = 'tomato')
      densinfos[[i]]$cond.dens[[j]] <- dens
      densinfos[[i]]$dnorm.dens[[j]] <- dn
      densinfos[[i]]$ylims[[j]] <- ylims
      densinfos[[i]]$xs[[j]] <- xs
    }
  }
  par(mfrow=c(1,1))
  dev.off()
  cat(sprintf("%s panel plot done.\n", toupper(family)))
  return(densinfos)
  on.exit(gc())
}

G <- plotFamily('gaussian')
P <- plotFamily('poisson')
B <- plotFamily('binomial')

source('zzz/theme_csda.R')
allPlots <- function(G, P, B, save.dir = './paper-sims/MVNjustification/output/'){
  # Gaussian ----------------------------------------------------------------
  Gr5b0 <- data.frame(condx = G$r5$cond.dens$b0$x, condy = G$r5$cond.dens$b0$y,
                      normy = G$r5$dnorm.dens$b0,
                      family = 'Gaussian', var = "b[0]", r = 5)
  Gr5b1 <- data.frame(condx = G$r5$cond.dens$b1$x, condy = G$r5$cond.dens$b1$y,
                      normy = G$r5$dnorm.dens$b1,
                      family = 'Gaussian', var = "b[1]", r = 5)
  Gr10b0 <- data.frame(condx = G$r10$cond.dens$b0$x, condy = G$r10$cond.dens$b0$y,
                       normy = G$r10$dnorm.dens$b0,
                       family = 'Gaussian', var = "b[0]", r = 10)
  Gr10b1 <- data.frame(condx = G$r10$cond.dens$b1$x, condy = G$r10$cond.dens$b1$y,
                       normy = G$r10$dnorm.dens$b1,
                       family = 'Gaussian', var = "b[1]", r = 10)
  Gr15b0 <- data.frame(condx = G$r15$cond.dens$b0$x, condy = G$r15$cond.dens$b0$y,
                       normy = G$r15$dnorm.dens$b0,
                       family = 'Gaussian', var = "b[0]", r = 15)
  Gr15b1 <- data.frame(condx = G$r15$cond.dens$b1$x, condy = G$r15$cond.dens$b1$y,
                       normy = G$r15$dnorm.dens$b1,
                       family = 'Gaussian', var = "b[1]", r = 15)
  Gdf <- rbind(Gr5b0, Gr5b1, Gr10b0, Gr10b1, Gr15b0, Gr15b1)
  # Poisson ----------------------------------------------------------------
  Pr5b0 <- data.frame(condx = P$r5$cond.dens$b0$x, condy = P$r5$cond.dens$b0$y,
                      normy = P$r5$dnorm.dens$b0,
                      family = 'Poisson', var = "b[0]", r = 5)
  Pr5b1 <- data.frame(condx = P$r5$cond.dens$b1$x, condy = P$r5$cond.dens$b1$y,
                      normy = P$r5$dnorm.dens$b1,
                      family = 'Poisson', var = "b[1]", r = 5)
  Pr10b0 <- data.frame(condx = P$r10$cond.dens$b0$x, condy = P$r10$cond.dens$b0$y,
                       normy = P$r10$dnorm.dens$b0,
                       family = 'Poisson', var = "b[0]", r = 10)
  Pr10b1 <- data.frame(condx = P$r10$cond.dens$b1$x, condy = P$r10$cond.dens$b1$y,
                       normy = P$r10$dnorm.dens$b1,
                       family = 'Poisson', var = "b[1]", r = 10)
  Pr15b0 <- data.frame(condx = P$r15$cond.dens$b0$x, condy = P$r15$cond.dens$b0$y,
                       normy = P$r15$dnorm.dens$b0,
                       family = 'Poisson', var = "b[0]", r = 15)
  Pr15b1 <- data.frame(condx = P$r15$cond.dens$b1$x, condy = P$r15$cond.dens$b1$y,
                       normy = P$r15$dnorm.dens$b1,
                       family = 'Poisson', var = "b[1]", r = 15)
  Pdf <- rbind(Pr5b0, Pr5b1, Pr10b0, Pr10b1, Pr15b0, Pr15b1)
  # Binomial ----------------------------------------------------------------
  # Remember, intercept only!
  Br5b0 <- data.frame(condx = B$r5$cond.dens$b0$x, condy = B$r5$cond.dens$b0$y,
                      normy = B$r5$dnorm.dens$b0,
                      family = 'Binomial', var = "b[0]", r = 5)
  Br10b0 <- data.frame(condx = B$r10$cond.dens$b0$x, condy = B$r10$cond.dens$b0$y,
                       normy = B$r10$dnorm.dens$b0,
                       family = 'Binomial', var = "b[0]", r = 10)
  Br15b0 <- data.frame(condx = B$r15$cond.dens$b0$x, condy = B$r15$cond.dens$b0$y,
                       normy = B$r15$dnorm.dens$b0,
                       family = 'Binomial', var = "b[0]", r = 15)
  Bdf <- rbind(Br5b0, Br10b0, Br15b0)
  # Make 3 ggplots and then combine ---------------------------------------
  fn <- paste0(save.dir, 'all.png')
  Gplot <- ggplot(Gdf, aes(x = condx, y = condy, colour = as.factor(r), group = as.factor(r))) + 
    geom_line() + 
    geom_line(aes(y = normy, col = as.factor(r)), lty = 5) + 
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
    theme_csda() + 
    labs(y = bquote(underline('Gaussian')),
         x = '', colour = expression(r~"=")) + theme(legend.position = 'none')
  Pplot <- ggplot(Pdf, aes(x = condx, y = condy, colour = as.factor(r), group = as.factor(r))) + 
    geom_line() + 
    geom_line(aes(y = normy, col = as.factor(r)), lty = 5) + 
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
    theme_csda() + 
    labs(y = bquote(underline('Poisson')),
         x = '', colour = expression(r~"=")) + theme(legend.position = 'none')
  Bplot <- ggplot(Bdf, aes(x = condx, y = condy, colour = as.factor(r), group = as.factor(r))) + 
    geom_line() + 
    geom_line(aes(y = normy, col = as.factor(r)), lty = 5) + 
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
    theme_csda() + 
    labs(y = bquote(underline('Binomial')),
         x = '', colour = expression(m[i]~"=")) + 
    theme(legend.position = 'none')
  png(fn,width = 140, height = 120, units = 'mm', pointsize = 9, res = 2000)
    gridExtra::grid.arrange(Gplot, Pplot, Bplot, nrow = 3)
  dev.off()
  cat("All-in-one plot done\n")
  on.exit(gc())
}
allPlots(G,P,B)
