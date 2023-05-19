source('zzz/theme_csda.R')
library(ggplot2)
library(dplyr)
theta <- c(-1, 0.0)
.sim <- function(family, n = 100){
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
  a <- simData(n = n, ntms = 15, theta = theta,
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
Sample <- function(data, btrue, family, ids, D){
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
  surv <- parseCoxph(Surv(survtime, status) ~ bin, data)
  l0 <- exp(theta[1] + theta[2] * fts)
  if(family == "binomial"){
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1|id)), parseFormula), l0)
    b.inds <- list(0)
    gamma.rep <- 0.5
  }else if(family == "gaussian"){
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
    b.inds <- list(0:1)
    gamma.rep <- c(0.5, 0.5)
  }else{
    sv <- surv.mod(surv, lapply(list(Y.1~time + cont + bin + (1 + time|id)), parseFormula), l0)  
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
  
  cli::cli_progress_bar(name = "Sampling...", total = length(ids))
  norm.dens <- cond.dens <- MVN.dens <- vector("list", length(ids))
  q <- ncol(D)
  Acc <- numeric(length(ids))
  for(a in 1:length(ids)){
    cond.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    norm.dens[[a]] <- setNames(vector("list", q), paste0("b", 0:(q-1)))
    # Walks
    W <- metropolis(b[[a]], Omega, Y[[a]], X[[a]], Z[[a]], 
                    list(family), Delta[[a]], S[[a]], Fi[[a]], l0i[[a]], 
                    SS[[a]], Fu[[a]], l0u[[a]], gamma.rep,
                    list(0:3), b.inds, 1L, length(b.inds[[1]]), 1000, 10000, Sigma[[a]], tune)
    Acc[a] <- W$AcceptanceRate
    for(j in 1:q){
      cond.dens[[a]][[j]] <- density(t(W$walks)[,j])
      norm.dens[[a]][[j]] <- dnorm(cond.dens[[a]][[j]]$x, mean = b.hat[[a]][j], 
                                   sd = sqrt(Sigma[[a]][j,j]))
    }
    MVN.dens[[a]] <- mvtnorm::dmvnorm(t(W$walks), mean = b.hat[[a]], sigma = Sigma[[a]])
    rm(W)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  # Make a data.frame to export.
  mi <- sapply(Z, function(x) nrow(x[[1]]))
  
  dfs <- lapply(1:length(ids), function(i){
    id <- i
    m <- mi[i]
    b0 <- data.frame(condx = cond.dens[[i]]$b0$x, condy = cond.dens[[i]]$b0$y,
                        normy = norm.dens[[i]]$b0,
                        family = family, var = "b[0]")
    if(q > 1){
      b1 <- data.frame(condx = cond.dens[[i]]$b1$x, condy = cond.dens[[i]]$b1$y,
                       normy = norm.dens[[i]]$b1,
                       family = family, var = "b[1]")
      i.dens <- rbind(b0,b1)
    }else{
      i.dens <- b0
    }
    i.dens$id <- i; i.dens$m <- mi[i]
    i.dens
  })
  
  out <- do.call(rbind, dfs) # About 4MB at 10,000 iterations on 100 subjects.
  out$ApproxBias <- out$normy - out$condy
  list(df = out, MVN.dens = MVN.dens, mi = unname(mi),
       family = family,
       Acc = Acc)
}

getOUT <- function(n, family){ # Wrapper for simulation + Sampling
  d <- .sim(family, n)
  btrue <- d[[2]]; D <- d[[3]]; data <- d[[1]]
  OUT <- Sample(data, btrue, family, 1:n, D)
}

plotOut <- function(OUT, save.dir = './paper-sims/MVNjustification/output/'){
  fn <- paste0(save.dir, OUT$family, "_RE_mi_breakdown.png")
  pp <- quantile(OUT$Acc, probs = c(0.025, 0.500, 0.975))
  cat(sprintf("Median [2.5%%, 97.5%%] acceptance rate for %s: %.2f [%.2f, %.2f]\n\n",
              OUT$family, pp[2], pp[1], pp[3]))
  df <- OUT$df
  
  # Sort into groups / make labels
  mi <- OUT$mi
  grp <- ifelse(mi <= 5, "[1,5]", ifelse(mi <= 10, "[6,10]", ifelse(mi <= 15,"[11,15]", "zzz")))
  grp <- factor(grp, 
                levels = c("[1,5]", "[6,10]", "[11,15]"),
                labels = c('"["*1*","*5*"]"', '"["*6*","*10*"]"', '"["*11*","*15*"]"'))
  
  mi.tab <- table(mi)
  mi.lab <- paste0(names(mi.tab[mi]), "\n(n=", mi.tab[mi],')')
  grp.tab <- table(grp)
  grp.lab <- paste0("atop(m[i] %in% ", names(grp.tab[grp]), ',"("*n==', grp.tab[grp],'*")")')
  
  mini <- data.frame(m = mi, mi.grp = grp, mi.lab = mi.lab, mi.grp.lab = grp.lab)
  mini <- mini[!duplicated(mini),]
  mini <- mini[order(mini$m),]
  mini$mi.grp <- forcats::fct_inorder(mini$mi.grp)
  mini$mi.lab <- forcats::fct_inorder(mini$mi.lab)
  mini$mi.grp.lab <- forcats::fct_inorder(mini$mi.grp.lab)
  
  df <- left_join(df, mini, 'm')
  
  # Plot (lots of) densities
  if(OUT$family != "binomial"){
    b0plot <- df %>% filter(var=='b[0]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_grid(~mi.grp.lab, scales = 'free', labeller = label_parsed)+
      labs(y = bquote(b[0]),x='') + 
      theme_csda()
    b1plot <- df %>% filter(var=='b[1]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_grid(~mi.grp.lab, scales = 'free', labeller = label_parsed)+
      labs(y = bquote(b[1]), x = '') + 
      theme_csda()+
      theme(
        strip.text = element_blank() # Remove top text as they're stacked.
      )
    png(fn, width = 190, height = 120, units = 'mm', pointsize = 9, res = 1000)
    gridExtra::grid.arrange(b0plot,b1plot, nrow=2,ncol=1)
    dev.off()
  }else{
    b0plot <- df %>% filter(var=='b[0]') %>% 
      ggplot(aes(x = condx, y = condy, group = id)) + 
      geom_line(lwd = .5) + 
      geom_line(aes(y = normy), lwd = .5, col = 'tomato', lty = 5) +
      facet_grid(~mi.grp.lab, scales = 'free', labeller = label_parsed)+
      labs(y = bquote(b[0]), x = '') + 
      theme_csda()
    png(fn, width = 190, height = 120, units = 'mm', pointsize = 9, res = 1000)
    print(b0plot)
    dev.off()
  }
  
  cat("Done for", OUT$family, ".\n")
}

getPlot <- function(n, family){
  OUT <- getOUT(n, family)
  plotOut(OUT)
  rm(OUT)
  on.exit(gc())
}

theta <- c(-1, 0.0)
getPlot(30, "gaussian")
getPlot(30, "poisson")
theta <- c(-1, 0.1)
getPlot(30, "binomial")
