# Quadratic vs linear time specification for serbilir ---------------------
rm(list=ls())
data(PBC)
PBC$serBilir <- log(PBC$serBilir)
long.formulas1 <- list(serBilir ~ drug * time + (1 + time|id))
long.formulas2 <- list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
surv.formula <- Surv(survtime, status) ~ drug
family <- list('gaussian')
fit <- joint(long.formulas1, surv.formula, PBC, family, control = list(verbose=T))
fit2 <- joint(long.formulas2, surv.formula, PBC, family, control = list(verbose=T))

logLik(fit)
logLik(fit2) # Appears to fit data better.

AIC(fit)
AIC(fit2) # Here too.

# What about AUCs?
fits <- list(fit, fit2)
