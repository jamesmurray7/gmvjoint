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

test <- PBC[PBC$id == 81 & PBC$time <= 6.9,]
ft <- my.fit$hazard[,1]
u <- ft[ft <= 15 & ft > 6.9]
jmlsurv <- joineRML::dynSurv(jML.fit, test, u = u, type = 'simulated', M = 200)
mine.normal <- dynPred(data = PBC, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 'normal')
mine.normal2 <- dynPred(data = PBC, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 'normal', scale = 2)
mine.t <- dynPred(data = PBC, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 't', df = 4, scale = 2)

plot(jmlsurv)
plot(mine.normal)




# Example from Riz book (p.176)
joineRML::dynSurv(jML.fit, PBC[PBC$id == 2,], type = 'simulated')
dynPred(PBC, 2, my.fit)
