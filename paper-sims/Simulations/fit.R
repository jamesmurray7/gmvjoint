long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin

fit <- function(d, ...){ # d is a single item from e.g. sim.sets$`...` .
  this.fit <- tryCatch(
    joint(long.formulas, surv.formula, d, 
          family = list("gaussian", "poisson", "binomial"),
          ...),
    error = function(e) NULL
  )
  this.fit
}

.pb <- function() utils::txtProgressBar(max = 100, style = 3)

fitwrap <- function(x){ # x is a list of data sets e.g. sim.sets$`mi = 5, failure = low`.
  out <- vector('list', 100)
  pb <- .pb()
  for(j in 1:100){
    out[[j]] <- fit(x[[j]], control = list(return.dmats = FALSE))
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  cat("\n")
  out
}

cat(names(sim.sets),sep='\n')
fits <- setNames(lapply(sim.sets, fitwrap), names(sim.sets))
save(fits, file='~/Downloads/fits.RData')


# Five-variate ------------------------------------------------------------

long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1 + time|id),
  Y.4 ~ time + cont + bin + (1 + time|id),
  Y.5 ~ time + cont + bin + (1|id)
)

surv.formula <- Surv(survtime, status) ~ bin

fit.five <- function(d, ...){ # d is a single item from e.g. sim.sets$`...` .
  this.fit <- tryCatch(
    joint(long.formulas, surv.formula, d, 
          family = list("gaussian", "gaussian", "poisson", "poisson", "binomial"),
          ...),
    error = function(e) NULL
  )
  this.fit
}

pb <- utils::txtProgressBar(max = 500, style = 3)
fits <- vector("list", 500)
for(j in 1:500){
  d <- sim.sets[[j]]
  out <- fit.five(d, control = list(return.dmats = FALSE))
  fits[[j]] <- out
  utils::setTxtProgressBar(pb, j)
}
save(fits, file = "~/Downloads/fits5variate.RData")