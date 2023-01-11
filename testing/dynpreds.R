# Testing probabilities I get are in keeping with other packages ----------
data(PBC)
PBC$serBilir <- log(PBC$serBilir)

#' Univariate
long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ drug
family <- list('gaussian')
my.fit <- joint(long.formulas, surv.formula, PBC, family, control = list(verbose=T))

library(joineRML)
jML.fit <- mjoint(
  formLongFixed = list('1' = serBilir ~ drug * time),
  formLongRandom = list('1' = ~ time | id),
  formSurv = Surv(survtime, status) ~ drug,
  data = PBC, timeVar = 'time', control = list(
    type = 'sobol', convCrit = 'rel', tol2 = 1e-1, tol.em = 5e-3
  ), verbose = F
)

test <- PBC[PBC$id == 81,]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = NULL, type = 'simulated', M = 100)
mine.normal <- dynPred(data = PBC, id = 81, fit = my.fit, u = NULL, nsim = 100, b.density = 'normal')
mine.normal2 <- dynPred(data = PBC, id = 81, fit = my.fit, u = NULL, nsim = 100, b.density = 'normal', scale = 0.33)
mine.t <- dynPred(data = PBC, id = 81, fit = my.fit, u = NULL, nsim = 100, b.density = 't', df = 4, scale = 2)

test <- PBC[PBC$id == 2,]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = NULL, type = 'simulated', M = 200)
mine.normal <- dynPred(data = PBC, id = 2, fit = my.fit, u = NULL, nsim = 200, b.density = 'normal', scale = 0.33)
mine.t <- dynPred(data = PBC, id = 2, fit = my.fit, u = NULL, nsim = 200, b.density = 't', df = 4, scale = 2)

plot(jmlsurv)
plot(mine.normal)

# Multivariate? -----------------------------------------------------------
long.formulas <- list(
  serBilir ~ drug * time + I(time^2) + (1 + time + I(time^2)|id),
  albumin ~ drug * time + (1 + time|id)
)

my.fit <- joint(long.formulas, surv.formula, data = PBC, family = list("gaussian", "gaussian"),
                control = list(verbose = T))

jML.fit <- mjoint(
  formLongFixed = list('1' = serBilir ~ drug * time + I(time^2),
                       '2' = albumin ~ drug * time),
  formLongRandom = list('1' = ~ time + I(time^2)|id, '2' = ~ time | id),
  formSurv = Surv(survtime, status) ~ drug,
  data = PBC, timeVar = 'time', control = list(
    type = 'sobol', convCrit = 'rel', tol2 = 1e-1, tol.em = 5e-3
  ), verbose = F
)

test <- PBC[PBC$id == 81,]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = NULL, type = 'simulated', M = 100, scale = 5)
mine.normal <- dynPred(data = PBC, id = 81, fit = my.fit, u = NULL, nsim = 100, b.density = 'normal', scale = 0.5)

plot(jmlsurv)
plot(mine.normal)
