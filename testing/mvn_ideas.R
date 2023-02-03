# Fit Biv. ----------------------------------------------------------------
data(PBC)
PBC <- na.omit(PBC[,c('id', 'albumin', 'time', 'survtime', 'status', 'platelets', 'drug' )])
PBC
long.formulas <- list(albumin~time*drug+(1+time|id),platelets~time*drug+(1+time|id))
surv.formula <- Surv(survtime, status)~drug
devtools::load_all(".")
fit <- joint(long.formulas, surv.formula, PBC, family = list("gaussian", "poisson"))


# ggplots -----------------------------------------------------------------
library(ggplot2)
theme_set(theme_light())

# 2d density plots --------------------------------------------------------
ggplot(data = (as.data.frame(fit$REs[,1:2])), 
       aes(x = `albumin_(Intercept)`, y = `albumin_time`)) + 
  geom_density_2d(data = as.data.frame(MASS::mvrnorm(10000, mu = c(0,0), 
                                                     Sigma = fit$coeffs$D[1:2,1:2])),
                  aes(`albumin_(Intercept)`, y = `albumin_time`),
                  colour = 'red', alpha = .5) + 
       geom_point(pch = 20, alpha = .5) + geom_density_2d()
  

# stat_ellipse ------------------------------------------------------------
ggplot(data = (as.data.frame(fit$REs[,1:2])), 
       aes(x = `albumin_(Intercept)`, y = `albumin_time`)) + 
  stat_ellipse(data = as.data.frame(MASS::mvrnorm(10000, mu = c(0,0), 
                                                     Sigma = fit$coeffs$D[1:2,1:2])),
                  aes(`albumin_(Intercept)`, y = `albumin_time`),
                  colour = 'red', alpha = .5, type = 'norm', lty = 2) + 
  stat_ellipse(aes(`albumin_(Intercept)`, y = `albumin_time`)) +
  geom_point(pch = 20, alpha = .5) 


# Singular density with theoretical curve? --------------------------------
par(mfrow = c(2,2))
for(q in 1:4){
  dens <- density(fit$REs[,q])
  plot(dens$x, dens$y, main = bquote("b"[.(q)]), type = 'l', xlab = 'Estimate',
       ylab = "density")
  curve(dnorm(x, mean = 0, sd = sqrt(fit$coeffs$D[q,q])), from = min(dens$x), to = max(dens$x),
        add = TRUE, col = 'red', lty = 'dashed')
}
par(mfrow = c(1,1))


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# How do appx. look with splines?
# Fit Spline? -------------------------------------------------------------
library(splines)
fit.spl <- joint(list(platelets ~ ns(time, 3) + (1 + ns(time, 3)|id)),
                 surv.formula, PBC, family = list('poisson'))

par(mfrow = c(2,2))
for(q in 1:4){
  dens <- density(fit.spl$REs[,q])
  plot(dens$x, dens$y, main = bquote("b"[.(q)]), type = 'l', xlab = 'Estimate',
       ylab = "density")
  curve(dnorm(x, mean = 0, sd = sqrt(fit.spl$coeffs$D[q,q])), from = min(dens$x), to = max(dens$x),
        add = TRUE, col = 'red', lty = 'dashed')
}
par(mfrow = c(1,1))



# What about a simulation? ------------------------------------------------
data <- simData(D = as.matrix(do.call(Matrix::bdiag, 
                            replicate(2, matrix(c(.25, 0.015, 0.015, .05), 2, 2), simplify = F))),
                theta = c(-2, .15), ntms = 2)$data

fit.sim <- joint(
  list(Y.1 ~ time + cont + bin + (1 + time|id), Y.2 ~ time + cont + bin + (1 + time|id)),
  Surv(survtime, status) ~ cont + bin,
  data,
  family = list("gaussian", "gaussian")
)

hist(fit.sim$dmats$surv$Tis[unlist(fit.sim$dmats$ph$Delta)==1L])

par(mfrow = c(2,2))
for(q in 1:4){
  dens <- density(fit.sim$REs[,q])
  plot(dens$x, dens$y, main = bquote("b"[.(q)]), type = 'l', xlab = 'Estimate',
       ylab = "density")
  curve(dnorm(x, mean = 0, sd = sqrt(fit.sim$coeffs$D[q,q])), from = min(dens$x), to = max(dens$x),
        add = TRUE, col = 'red', lty = 'dashed')
}
par(mfrow = c(1,1))


ggplot(data = (as.data.frame(fit.sim$REs[,1:2])), 
       aes(x = `Y.1_(Intercept)`, y = `Y.1_time`)) + 
  geom_density_2d(data = as.data.frame(MASS::mvrnorm(10000, mu = c(0,0), 
                                                     Sigma = fit.sim$coeffs$D[1:2,1:2])),
                  aes(`Y.1_(Intercept)`, y = `Y.1_time`),
                  colour = 'red', alpha = .5, lineend = 'round', linejoin = 'round') + 
  geom_point(pch = 20, alpha = .5) + geom_density_2d()

ggplot(data = (as.data.frame(fit.sim$REs[,1:2])), 
       aes(x = `Y.1_(Intercept)`, y = `Y.1_time`)) + 
  stat_ellipse(data = as.data.frame(MASS::mvrnorm(10000, mu = c(0,0), 
                                                     Sigma = fit.sim$coeffs$D[1:2,1:2])),
                  aes(`Y.1_(Intercept)`, y = `Y.1_time`),
                colour = 'red', alpha = .5, type = 'norm', lty = 2) +
  geom_point(pch = 20, alpha = .5) + stat_ellipse(type='norm') + 
  labs(y = bquote(b[1]), xlab = bquote(b[0]))
