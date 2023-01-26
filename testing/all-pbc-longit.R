library(splines)
srv <- Surv(survtime, status) ~ drug
# serBilir ----------------------------------------------------------------
lin <- list(serBilir ~ drug * time + (1 + time|id))
qua <- list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
spl <- list(serBilir ~ drug * (ns(time, 3)) + (1 + ns(time, 3)|id))


fit.lin <- joint(lin, srv, PBC, list("gaussian"))
fit.qua <- joint(qua, srv, PBC, list("gaussian"))
fit.spl <- joint(spl, srv, PBC, list("gaussian"))

anova(fit.lin, fit.qua) # p < 0.001
anova(fit.qua, fit.spl) # p = 0.647; 


# Prot --------------------------------------------------------------------
PBC$prot <- (0.1 * PBC$prothrombin)^(-4)
any(is.na(PBC$prot))

lin <- list(prot ~ drug * time + (1 + time|id))
qua <- list(prot ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
spl <- list(prot ~ drug * (ns(time, 3)) + (1 + ns(time, 3)|id))

fit.lin <- joint(lin, srv, PBC, list("gaussian"))
fit.qua <- joint(qua, srv, PBC, list("gaussian"))
fit.spl <- joint(spl, srv, PBC, list("gaussian"))

anova(fit.lin, fit.qua) # p < 0.001
anova(fit.qua, fit.spl) # p < 0.001
plot(resid(fit.spl))


# Albumin -----------------------------------------------------------------
any(is.na(PBC$albumin))

lin <- list(albumin ~ drug * time + (1 + time|id))
qua <- list(albumin ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
spl <- list(albumin ~ drug * (ns(time, 3)) + (1 + ns(time, 3)|id))

fit.lin <- joint(lin, srv, PBC, list("gaussian"))
fit.qua <- joint(qua, srv, PBC, list("gaussian"))
fit.spl <- joint(spl, srv, PBC, list("gaussian"))

anova(fit.lin, fit.qua) # p < 0.001
anova(fit.lin, fit.spl)
anova(fit.qua, fit.spl) # p < 0.001
plot(resid(fit.spl))    # Spline supposedly best? Think this might be due to glmmTMB poorly fitting quad.

glmmTMB(qua[[1]], PBC, REML = F)
plot(resid(fit.spl, type = 'pearson'))


# Platelets ---------------------------------------------------------------
data <- PBC[!is.na(PBC$platelets), ]

lin <- list(platelets ~ drug * time + (1 + time|id))
qua <- list(platelets ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
spl <- list(platelets ~ drug * (ns(time, 3)) + (1 + ns(time, 3)|id))

fit.lin <- joint(lin, srv, data, list("poisson"))
fit.qua <- joint(qua, srv, data, list("poisson"))
fit.spl <- joint(spl, srv, data, list("poisson"))

fit.lin; fit.qua; fit.spl

anova(fit.lin, fit.qua)
anova(fit.qua, fit.spl) # spline.





