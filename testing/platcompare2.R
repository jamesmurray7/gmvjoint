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
AUCsG <- bootAUC(fitG, PBC, Tstart = 8, delta = 2, nboot = 100)
AUCsP <- bootAUC(fitP, PBC, Tstart = 8, delta = 2, nboot = 100)
AUCsGP <- bootAUC(fitGP, PBC, Tstart = 8, delta = 2, nboot = 100)

AUCsG
AUCsP
AUCsGP

boxplot(AUCsG$AUCs, AUCsP$AUCs, AUCsGP$AUCs, xaxt = 'n')
axis(1, at = 1:3, labels = c('Gaussian', 'Poisson', 'Generalised Poisson'), tick = F)

# 100 first order estimates in window (4, 7]
AUCsG2 <- bootAUC(fitG, PBC, Tstart = 4, delta = 3, nboot = 100)
AUCsP2 <- bootAUC(fitP, PBC, Tstart = 4, delta = 3, nboot = 100)
AUCsGP2 <- bootAUC(fitGP, PBC, Tstart = 4, delta = 3, nboot = 100)

AUCsG2
AUCsP2
AUCsGP2

boxplot(AUCsG$AUCs, AUCsP$AUCs, AUCsGP$AUCs, xaxt = 'n')
axis(1, at = 1:3, labels = c('Gaussian', 'Poisson', 'Generalised Poisson'), tick = F)