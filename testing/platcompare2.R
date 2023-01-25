rm(list=ls())
devtools::load_all()
data(PBC)
PBC <- na.omit(PBC[,c('id', 'time', 'survtime', 'status', 'drug', 'platelets')])
PBC$plat2 <- c(scale(PBC$platelets))
long.formulas <- list(platelets ~ drug * time + (1 + time|id))
surv.formula <- Surv(survtime, status) ~ drug
fitG <- joint(list(plat2 ~ drug * time + (1 + time|id)), surv.formula, PBC, list('gaussian'), control = list(verbose=T))
fitP <- joint(long.formulas, surv.formula, PBC, list('poisson'), control = list(verbose=T))
fitGP <- joint(long.formulas, surv.formula, PBC, list('genpois'), control = list(verbose=T))
# 100 first order estimates in window (8, 10]
AUCsG <- bootAUC(fitG, PBC, Tstart = 8, delta = 2, nboot = 5)
AUCsP <- bootAUC(fitP, PBC, Tstart = 8, delta = 2, nboot = 5)
AUCsGP <- bootAUC(fitGP, PBC, Tstart = 8, delta = 2, nboot = 5)

AUCsG
AUCsP
AUCsGP
