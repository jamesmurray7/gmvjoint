# Simulate data -----------------------------------------------------------
N <- 100
datas <- replicate(N, 
                   simData(D = diag(c(.25, .05)), beta = t(c(2, 0.33, -0.5, 0.25)),
                           theta = c(-2.9, 0.1),
                           gamma = 0.5, family = list('gaussian'))$dat,
                   simplify = F)
# save(datas, file = "/data/c0061461/fits_ghnodes/datas.RData")

ff <- function(dat, gh){
  dd <- dat
  con <- list(gh.nodes = gh,
              tol.rel = 1e-3,
              return.dmats = FALSE)
  fit <- joint(list(Y.1 ~ time + cont + bin + (1 + time|id)),
               Surv(survtime, status) ~ bin, data = dd,
               family = list('gaussian'),
               control = con)
  fit
}

library(cli)
nodes <- c(1:3,5,7,9,15)

fit.gh <- function(gh, datas){
  nm <- paste0("Simulation study, number of GH nodes: ", gh)
  fn <- paste0("/data/c0061461/fits_ghnodes/n_", gh, '.RData')
  cli_progress_bar(name = nm, total = N)
  fits <- vector("list", N)
  for(i in 1:N){
    this <- tryCatch(suppressMessages(ff(datas[[i]], gh)),
                     error = function(e) NULL)
    if(is.null(this))
      fits[[i]] <- NULL
    else
      fits[[i]] <- this
    cli_progress_update()
  }
  cli_progress_done()
  cli_alert_success(sprintf("\nSaving in %s\n", fn))
  out <- fits
  save(out, file = fn)
}

for(g in nodes) fit.gh(g, datas)


# Trivariate mix. ---------------------------------------------------------
rm(list=ls())
N <- 100
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

N <- 100
datas <- replicate(N,  # appx. 30%
                   simData(n = 250, ntms = 10, family = family, sigma = sigma, beta = beta,
                           D = D, gamma = gamma, zeta = zeta, theta = c(-2.9, .1),
                           random.formula = list(~time, ~time, ~1))$data,
                   simplify = F)
save(datas, file = "/data/c0061461/fits_ghnodes/datas.RData")

ff <- function(dat, gh){
  dd <- dat
  con <- list(gh.nodes = gh,
              tol.rel = 1e-2,
              return.dmats = FALSE)
  fit <- joint(list(Y.1 ~ time + cont + bin + (1 + time|id),
                    Y.2 ~ time + cont + bin + (1 + time|id),
                    Y.3 ~ time + cont + bin + (1|id)),
               Surv(survtime, status) ~ bin, data = dd,
               family = list('gaussian', 'poisson', 'binomial'),
               control = con)
  fit
}

nodes <- c(2:3,5,7,9,15)

fit.gh <- function(gh, datas){
  nm <- paste0("Trivariate simulation study, number of GH nodes: ", gh)
  fn <- paste0("/data/c0061461/fits_ghnodes/n_", gh, '.RData')
  cli_progress_bar(name = nm, total = N)
  fits <- vector("list", N)
  for(i in 1:N){
    this <- tryCatch(suppressMessages(ff(datas[[i]], gh)),
                     error = function(e) NULL)
    if(is.null(this))
      fits[[i]] <- NULL
    else
      fits[[i]] <- this
    cli_progress_update()
  }
  cli_progress_done()
  cli_alert_success(sprintf("\nSaving in %s\n", fn))
  out <- fits
  save(out, file = fn)
}

for(g in nodes) fit.gh(g, datas)


# Trivariate Mix, bigger covariance ---------------------------------------
rm(list=ls())
N <- 100
D <- diag(c(2, .4, 1.5, 0.33, 3.00))
D[1,3] <- D[3,1] <- D[1,5] <- D[5,1] <- D[3,5] <- D[5,3] <- 0.5
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

N <- 100
datas <- replicate(N,  # appx. 30%
                   simData(n = 250, ntms = 10, family = family, sigma = sigma, beta = beta,
                           D = D, gamma = gamma, zeta = zeta, theta = c(-1.9, .1),
                           random.formula = list(~time, ~time, ~1))$data,
                   simplify = F)
save(datas, file = "/data/c0061461/fits_ghnodes_largeD/datas.RData")

ff <- function(dat, gh){
  dd <- dat
  con <- list(gh.nodes = gh,
              tol.rel = 1e-2,
              return.dmats = FALSE)
  fit <- joint(list(Y.1 ~ time + cont + bin + (1 + time|id),
                    Y.2 ~ time + cont + bin + (1 + time|id),
                    Y.3 ~ time + cont + bin + (1|id)),
               Surv(survtime, status) ~ bin, data = dd,
               family = list('gaussian', 'poisson', 'binomial'),
               control = con)
  fit
}

# https://arxiv.org/pdf/2202.07864.pdf
stringer.gh <- function(data, what = 'median'){
  if(what == "median")
    rmin <- median(with(data, tapply(time, id, length)))
  else
    rmin <- quantile(with(data, tapply(time, id, length)), probs = .25)
  n <- length(unique(data$id))
  if(rmin == 1) rmin <- rmin + 1
  ceiling(1.5 * log(n, base = rmin) - 2)
}

nodes <- c(2, 3, 5, 7, 9, 15, 98, 99)

fit.gh <- function(gh, datas){
  nm <- paste0("Trivariate simulation study, number of GH nodes: ", gh)
  fn <- paste0("/data/c0061461/fits_ghnodes_largeD/n_", gh, '.RData')
  cli_progress_bar(name = nm, total = N)
  fits <- vector("list", N)
  for(i in 1:N){
    if(gh == 98) gh <- stringer.gh(datas[[i]], 'median')
    if(gh == 99) gh <- stringer.gh(datas[[i]], 'qqqq')
    this <- tryCatch(suppressMessages(ff(datas[[i]], gh)),
                     error = function(e) NULL)
    if(is.null(this))
      fits[[i]] <- NULL
    else
      fits[[i]] <- this
    cli_progress_update()
  }
  cli_progress_done()
  cli_alert_success(sprintf("\nSaving in %s\n", fn))
  out <- fits
  save(out, file = fn)
}

for(g in nodes) fit.gh(g, datas)

