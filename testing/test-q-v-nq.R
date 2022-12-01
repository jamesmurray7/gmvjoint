long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ bin
family <- list('gaussian','poisson','binomial')
data <- simData(beta =  rbind(c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50)),
                family = family, gamma = c(.5,-.5,.5), sigma = rep(.2,3),
                D = diag(c(.25, .09,.25,.09,.25,.09)))
data <- data$data
control <- list()


fitnq <- joint(long.formulas, surv.formula, data, family, control = list(verbose = T))
summary(fitnq)

fitq <- joint(long.formulas, surv.formula, data, family, control = list(verbose = T, beta.quad=T))

summary(fitq)
