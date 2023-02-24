long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin

fit <- function(d, ...){
  this.fit <- tryCatch(
    joint(long.formulas, surv.formula, d, 
          family = list("gaussian", "poisson", "binomial"),
          ...),
    error = function(e) NULL
  )
  this.fit
}
