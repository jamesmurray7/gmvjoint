# Marginal distn. of b|Y for longit binomial distn
test <- function(b, Y, X, Z, beta, D){
  eta <- X %*% beta + Z %*% b
  prod(dbinom(Y,1,plogis(eta)) * mvtnorm::dmvnorm(b, sigma = D))
}

library(MCMCpack)
invest.ids <- sample(1:250, 5)

# look at Binom only...
A <- .sim(250, 5,  'binomial', fup = 15, unif.times = F)
B <- .sim(250, 10, 'binomial', fup = 15, unif.times = F)
C <- .sim(250, 15, 'binomial', fup = 15, unif.times = F)

fun1R <- function(dat){
  X <- Y <- Z <- vector('list', 250)
  for(i in 1:250){
    X[[i]] <- model.matrix(~time+cont+bin, dat[dat$id==i,,drop=F])
    Z[[i]] <- model.matrix(~time, dat[dat$id==i,,drop=F])
    Y[[i]] <- dat[dat$id==i,'Y.1']
  }
  
  Yl <- sapply(Y, length)
  qsamp <- quantile(Yl); med <- unname(qsamp[3]); max <- qsamp[5]
  med.to.samp <- sample(which(Yl == med),5,replace=T); max.to.samp <- sample(which(Yl == max)[1:5],5,replace=T)
  outmed <- setNames(vector('list', 5), paste0('med: mi = ', med, '_',1:5))
  outmax <- setNames(vector('list', 5), paste0('max: mi = ', max, '_',1:5))
  pb <- utils::txtProgressBar(max = 5, style = 3)
  for(j in 1:5){
    outmed[[j]] <- MCMCmetrop1R(test, theta.init = c(0,0), Y = Y[[med.to.samp[j]]], X = X[[med.to.samp[j]]], 
                                Z = Z[[med.to.samp[j]]], 
                                beta = c(0, -1, 0.5, -0.5), D = diag(c(.5, .1)),
                                burnin = 50, mcmc = 5000, logfun = F)
    outmax[[j]] <- MCMCmetrop1R(test, theta.init = c(0,0), Y = Y[[max.to.samp[j]]], X = X[[max.to.samp[j]]], 
                                Z = Z[[max.to.samp[j]]], 
                                beta = c(0, -1, 0.5, -0.5), D = diag(c(.5, .1)),
                                burnin = 50, mcmc = 5000, logfun = F)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  list(Med = outmed, Max = outmax)
}

Asamp <- fun1R(A)
Bsamp <- fun1R(B)
Csamp <- fun1R(C)

# Plotting ----------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
source('zzz/theme_csda.R')

plotSamps <- function(samps){
  meds <- samps[[1]]
  maxs <- samps[[2]]
  mi.med <- unique(gsub('\\_\\d$','',stringr::str_extract(names(meds), '\\d?\\d\\_\\d')))
  mi.max <- unique(gsub('\\_\\d$','',stringr::str_extract(names(maxs), '\\d?\\d\\_\\d')))
  
  # Function to get densities
  dn <- function(l){
    lapply(lapply(l, function(x) apply(x, 2, density)),
           function(y) data.frame(b0_x = y$var1$x, b0_dens = y$var1$y, b1_x = y$var2$x, b1_dens = y$var2$y))
  }
  med.df <- dn(meds); for(i in 1:length(med.df)) med.df[[i]]$id <- i;
  med.df <- do.call(rbind, med.df); med.df$m <- mi.med
  row.names(med.df) <- NULL
  max.df <- dn(maxs); for(i in 1:length(max.df)) max.df[[i]]$id <- (i + 5); 
  max.df <- do.call(rbind, max.df); max.df$m <- mi.max
  row.names(max.df) <- NULL
  
  # head(max.df); head(med.df)
  df <- rbind(med.df, max.df)
  fin <- df %>% 
    pivot_longer(b0_x:b1_dens) %>% 
    mutate(intslp = ifelse(grepl('b0', name), 'b[0]', 'b[1]'),
           what = ifelse(grepl('x$', name), 'x', 'y')) %>% 
    select(-name) %>% 
    pivot_wider(values_from = value, names_from = what, values_fn = list) %>% 
    unnest(c(x, y))

  fin$pred <- ifelse(fin$intslp=='b[0]', dnorm(fin$x, mean = mean(c(fin[fin$intslp=='b[0]', 'x'])$x), sd = sqrt(.5)),
                     dnorm(fin$x, mean = mean(c(fin[fin$intslp=='b[1]', 'x'])$x), sd = sqrt(.1)))
  
  fin %>% 
    ggplot(aes(x = x, y = y, group = id, colour = m)) + 
    geom_line() + 
    geom_line(aes(y = pred), colour = 'red', lwd = 5) + 
    facet_wrap(~intslp, scales = 'free', labeller = label_parsed) + 
    theme_csda()
    
  
}
plotSamps(Csamp)

# Plotting entire sample --------------------------------------------------
long.formulas <- list(Y.1~time+cont+bin+(1+time|id))
formulas <- lapply(long.formulas, parseFormula)
test <- function(b, Y, X, Z, beta, D, S, SS, Fu, Fi, l0i, haz, Del, gam, zet,
                 survival = TRUE){
  eta <- X %*% beta + Z %*% b
  out <- prod(dbinom(Y, 1, plogis(eta))) * mvtnorm::dmvnorm(b, sigma = D) 
  if(survival) out <- out * exp(logfti(b, S, SS, Fi, Fu, l0i, haz, Del, gam, zet))
  out
}

# 70% failure rates
A <- .sim(100, 5,  'binomial', fup = 15, unif.times = F)
B <- .sim(100, 10, 'binomial', fup = 15, unif.times = F)
C <- .sim(100, 15, 'binomial', fup = 15, unif.times = F)
E <- .sim(100, 30, 'binomial', fup = 15, unif.times = F)
fun1R <- function(dat, survival){
  # Longit
  X <- Y <- Z <- vector('list', 100)
  for(i in 1:100){
    X[[i]] <- model.matrix(~time+cont+bin, dat[dat$id==i,,drop=F])
    Z[[i]] <- model.matrix(~time, dat[dat$id==i,,drop=F])
    Y[[i]] <- dat[dat$id==i,'Y.1']
  }
  # Survival part
  fts <- sort(unique(dat[dat$status==1,'survtime']))
  surv <- parseCoxph(Surv(survtime, status) ~ bin, dat)
  l0 <- exp(-2.5 + 0.1 * fts)
  sv <- surv.mod(surv, formulas, l0)
  out <- Map(function(X, Y, Z, S, SS, Fu, Fi, l0i, haz, Del){
    MCMCmetrop1R(test, theta.init = c(0,0), 
                 Y = Y, X = X, Z = Z, beta = c(0, -1, 0.5, -0.5), D = diag(c(.5, .1)),
                 S = S, SS = SS, Fu = Fu, Fi = Fi, l0i = l0i,
                 haz = haz, Del = Del, gam = rep(-0.5, 2), zet = 0.3, survival = survival,
                 burnin = 500, mcmc = 3500, logfun = F, tune = 2.5)#c(.33, .33))
  }, X = X, Y = Y, Z = Z, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0i = sv$l0i,
     haz = sv$l0u, Del = surv$Delta)
}

# Quickly plot ONE 
singleplot <- function(x){
  par(mfrow=c(1,2))
  for(i in 1:2){
    d <- density(x[,i])
    plot(d, xlab = '', main = bquote(b[.(i-1)]))
    curve(dnorm(x, 0, sqrt(diag(c(.5,.1))[i,i])), from = min(d$x), to = max(d$x), add = T, 
          col = 'red', lty = 3)
  }
  par(mfrow=c(1,1))
}

compareplot <- function(x,y){
  par(mfrow=c(1,2))
  for(i in 1:2){
    dx <- density(x[,i]); dy <- density(y[,i])
    plot(dx, xlab = '', main = bquote(b[.(i-1)]),
         xlim = c(min(pmin(dx$x, dy$x)),max(pmax(dx$x, dy$x))),
         ylim = c(min(pmin(dx$y, dy$y)),max(pmax(dx$y, dy$y))))
    lines(dy, col = 'steelblue')
    curve(dnorm(x, 0, sqrt(diag(c(.5,.1))[i,i])), from = min(d$x), to = max(d$x), add = T, 
          col = 'red', lty = 3)
  }
  # Put legend on top-right only!
  legend('topright', bty = 'n', lty = c(1,1,3), col = c('black', 'steelblue', 'red'),
         legend = c(expression(f*"("*b[i]*"|"*Y[i]*";"*Omega*")"),
                    expression(f*"("*b[i]*"|"*Y[i]*","*T[i]*","*Delta[i]*";"*Omega*")"),
                    "Theoretical"))
  par(mfrow=c(1,1))
}

# mi = 5 -=-=-=
Asampfull <- fun1R(A, TRUE)
Asamplong <- fun1R(A, FALSE)
singleplot(do.call(rbind, Asamplong))
singleplot(do.call(rbind, Asampfull))
compareplot(do.call(rbind, Asamplong),do.call(rbind, Asampfull))

# mi = 15 -=-=-=
Csampfull <- fun1R(C, T)
Csamplong <- fun1R(C, F)
compareplot(do.call(rbind, Csamplong),do.call(rbind, Csampfull))

# 40-50% failure rates
A <- .sim(100, 5,  'binomial')
C <- .sim(100, 15, 'binomial')
Asampfull <- fun1R(A, TRUE)
Asamplong <- fun1R(A, FALSE)
singleplot(do.call(rbind, Asamplong))
singleplot(do.call(rbind, Asampfull))
compareplot(do.call(rbind, Asamplong),do.call(rbind, Asampfull))

# mi = 15 -=-=-=
Csampfull <- fun1R(C, T)
Csamplong <- fun1R(C, F)
compareplot(do.call(rbind, Csamplong),do.call(rbind, Csampfull))

# 20%ish failure rates
A <- simData(n = 100, ntms = 5, beta = t(c(0, -1, 0.5, -0.5)),
             D <- diag(c(.5, .1)), gamma = -.5, zeta = c(0, .3),
             family = list('binomial'), fup = 5,
             sigma=c(.16), theta = c(-3.5, .1))$data
C <- simData(n = 100, ntms = 15, beta = t(c(0, -1, 0.5, -0.5)),
             D <- diag(c(.5, .1)), gamma = -.5, zeta = c(0, .3),
             family = list('binomial'), fup = 5,
             sigma=c(.16), theta = c(-3.5, .1))$data
Asampfull <- fun1R(A, TRUE)
Asamplong <- fun1R(A, FALSE)
singleplot(do.call(rbind, Asamplong))
singleplot(do.call(rbind, Asampfull))
compareplot(do.call(rbind, Asamplong),do.call(rbind, Asampfull))

# mi = 15 -=-=-=
Csampfull <- fun1R(C, T)
Csamplong <- fun1R(C, F)
compareplot(do.call(rbind, Csamplong),do.call(rbind, Csampfull))
