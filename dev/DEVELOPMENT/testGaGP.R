family <- list("Gamma", "negbin")
sigma <- list(log(2), 0.5)
data <- simData(family = family,
                fup = 1, ntms = 10, beta = do.call(rbind, replicate(2, c(2, -0.5, 0.15, 0.2), FALSE)),
                disp.formulas = list(~time, ~time),
                sigma = list(c(2, -0.25), c(-1, 0.3)),
                gamma = rep(2, 2))$data
control <- list()
disp.formulas <- list(~time, ~time)
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ cont + bin
fit <- joint(long.formulas, surv.formula, disp.formulas = list(~time, ~time),
             data, family = list("Gamma", "negbin"),
             control = list(verbose=T))
