library(MCMCpack)
long.formulas <- list(Y.1~time+cont+bin+(1+time|id))
formulas <- lapply(long.formulas, parseFormula)
test <- function(b, Y, X, Z, beta, sigma, D, S, SS, Fu, Fi, l0i, haz, Del, gam, zet,
                 survival = TRUE){
  eta <- X %*% beta + Z %*% b
  out <- prod(dnorm(Y, eta, sd = sqrt(sigma))) * mvtnorm::dmvnorm(b, sigma = D) 
  if(survival) out <- out * exp(logfti(b, S, SS, Fi, Fu, l0i, haz, Del, gam, zet))
  out
}

# 70% failure rates
A <- .sim(100, 5,  'gaussian', fup = 15, unif.times = F)
B <- .sim(100, 10, 'gaussian', fup = 15, unif.times = F)
C <- .sim(100, 15, 'gaussian', fup = 15, unif.times = F)
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
    a <- tryCatch(MCMCmetrop1R(test, theta.init = c(0,0), 
                 Y = Y, X = X, Z = Z, beta = c(2, -.1, 0.1, .2), sigma = .16, D = diag(c(.25, .05)),
                 S = S, SS = SS, Fu = Fu, Fi = Fi, l0i = l0i,
                 haz = haz, Del = Del, gam = rep(-0.5, 2), zet = 0.3, survival = survival,
                 burnin = 500, mcmc = 3500, logfun = F, tune = 2.5),
                 error = function(e) NA
    )
    a
  }, X = X, Y = Y, Z = Z, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, l0i = sv$l0i,
  haz = sv$l0u, Del = surv$Delta)
}

# Quickly plot ONE 
singleplot <- function(x){
  par(mfrow=c(1,2))
  for(i in 1:2){
    d <- density(x[,i])
    plot(d, xlab = '', main = bquote(b[.(i-1)]))
    curve(dnorm(x, 0, sqrt(diag(c(.25,.05))[i,i])), from = min(d$x), to = max(d$x), add = T, 
          col = 'red', lty = 3)
  }
  par(mfrow=c(1,1))
}

compareplot <- function(x,y){
  par(mfrow=c(1,2))
  for(i in 1:2){
    dx <- density(x[,i], na.rm=T); dy <- density(y[,i],na.rm=T)
    plot(dx, xlab = '', main = bquote(b[.(i-1)]),
         xlim = c(min(pmin(dx$x, dy$x)),max(pmax(dx$x, dy$x))),
         ylim = c(min(pmin(dx$y, dy$y)),max(pmax(dx$y, dy$y))))
    lines(dy, col = 'steelblue')
    curve(dnorm(x, 0, sqrt(diag(c(.25,.05))[i,i])), from = min(d$x), to = max(d$x), add = T, 
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
A <- .sim(100, 5,  'gaussian')
C <- .sim(100, 15, 'gaussian')
Asampfull <- fun1R(A, TRUE)
Asamplong <- fun1R(A, FALSE)
compareplot(do.call(rbind, Asamplong),do.call(rbind, Asampfull))

# mi = 15 -=-=-=
Csampfull <- fun1R(C, T)
Csamplong <- fun1R(C, F)
compareplot(do.call(rbind, Csamplong),do.call(rbind, Csampfull))

# 20%ish failure rates
A <- simData(n = 100, ntms = 5, beta = t(c(2, -.1, 0.1, .2)),
             D <- diag(c(.25, .05)), gamma = -.5, zeta = c(0, .3),
             family = list('gaussian'), fup = 5,
             sigma=c(.16), theta = c(-3.5, .1))$data
C <- simData(n = 100, ntms = 15, beta = t(c(2, -.1, 0.1, .2)),
             D <- diag(c(.25, .05)), gamma = -.5, zeta = c(0, .3),
             family = list('gaussian'), fup = 5,
             sigma=c(.16), theta = c(-3.5, .1))$data
Asampfull <- fun1R(A, TRUE)
Asamplong <- fun1R(A, FALSE)
compareplot(do.call(rbind, Asamplong),do.call(rbind, Asampfull))

# mi = 15 -=-=-=
Csampfull <- fun1R(C, T)
Csamplong <- fun1R(C, F)
compareplot(do.call(rbind, Csamplong),do.call(rbind, Csampfull))
