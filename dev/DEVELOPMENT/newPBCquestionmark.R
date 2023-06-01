rm(list=ls())
source(".Rprofile")
library(splines)
library(xtable)
data(PBC, package = 'gmvjoint')

# Some issues with Ascites -- Cov matrix seems fairly volatile.
unlist(with(PBC, tapply(ascites, id, function(x) table(unique(x)))))
table(with(PBC, tapply(ascites, id, function(x){
  sum(x==1)/length(x)
})))
sum(PBC$ascites==1,na.rm=T)/nrow(gmvjoint::PBC)
sum(PBC$hepatomegaly==1)/nrow(PBC)
table(with(PBC, tapply(hepatomegaly, id, function(x){
  sum(x==1)/length(x)
})))
sum(PBC$spiders==1)/nrow(PBC)
table(with(PBC, tapply(spiders, id, function(x){
  sum(x==1)/length(x)
})))

PBC <- na.omit(PBC[,c("id", "survtime", "status", "drug", "sex", "age", "time",
                      "hepatomegaly", "spiders", "serBilir",
                      "albumin", "alkaline", "SGOT", "platelets", "prothrombin")])
PBC$serBilir <- log(PBC$serBilir)
# PBC$prothrombin <- (PBC$prothrombin * .1)^ (-4)
PBC$AST <- log(PBC$SGOT)

Gaussian.long.formulas <- list(
  serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
  albumin ~ drug * time + (1 + time|id),
  # prothrombin ~ drug * ns(time, knots = c(1, 4)) + (1 + ns(time,knots = c(1, 4))|id),
  prothrombin ~ drug * time + (1 + time|id),
  # AST ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id)
  AST~drug * time + (1 + time|id)
)

Poisson.long.formulas <- list(
  # platelets ~ drug * time  + (1 + time|id),
  platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1, 4))|id),
  alkaline ~ drug * time  + (1 + time|id)
)

Binomial.long.formulas <- list(
  hepatomegaly ~ drug * time  + (1|id),
  spiders ~ drug * time  + (1|id)
)

surv.formula <- Surv(survtime, status) ~ drug
control <- list(verbose=T)

# 4 variate Gaussian
Gaussians <- joint(Gaussian.long.formulas, surv.formula, PBC, 
                   list("gaussian", "gaussian", "gaussian", "gaussian"))

summary(Gaussians)
xtable(Gaussians, vcov = TRUE, max.row = 14, size = "footnotesize", booktabs = FALSE)
# Take forward: serBilir, albumin and AST

Poissons <- joint(Poisson.long.formulas,
                  surv.formula, PBC, list("poisson", "poisson"))

xtable(Poissons, vcov = TRUE, size = 'footnotesize', max.row = 10, booktabs = FALSE)

Binomials <- joint(Binomial.long.formulas, surv.formula,
                   PBC, list("binomial", "binomial"))
xtable(Binomials, vcov = TRUE, size = 'footnotesize', max.row = 7, booktabs = FALSE)

# Take Gaussian: {serBilir, albumin, prothrombin}, 
# Poisson: {platelets} and Binomial {hepatomegaly} forward...
reduced.long <- list(
  serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
  albumin ~ drug * time + (1 + time|id),
  AST ~ drug * time  + (1 + time|id),
  # prothrombin ~ drug * ns(time, 3) + (1 + ns(time, 3)|id),
  # prothrombin ~ ti
  # platelets ~ drug * time  + (1 + time|id),
  platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1, 4))|id),
  hepatomegaly ~ drug * time  + (1|id)
)

reduced.model <- joint(reduced.long,
                       surv.formula, PBC, list("gaussian", "gaussian", "gaussian",
                                               "poisson", "binomial"))

summary(reduced.model)
xtable(reduced.model, vcov = TRUE, size = 'footnotesize',
       booktabs = FALSE)


# JMbayes2 fit ------------------------------------------------------------
pt <- proc.time()[3]                                     # Start timing
lsb <- lme(fixed = serBilir ~ drug * (time + I(time^2)), # log(serBilir)
           random = ~ time + I(time^2)|id,
           data = PBC)
alb <- lme(fixed = albumin ~ drug * time,                # albumin
           random = ~ time|id,
           data = PBC)
ast <- lme(fixed = AST ~ drug * time,   # AST
           random = ~ time | id,
           data = PBC)
plt <- mixed_model(fixed = platelets ~ drug * ns(time, knots = c(1,4)),      # platetelet count
                   random = ~1 + ns(time, knots = c(1, 4))|id,               # mixed_model struggles to fit for absolutely no reason, 
                   data = PBC, family = poisson(),                           # so start at MLEs
                   initial_values = 
                     list(betas = Poissons$coeffs$beta[1:8],
                          D = Poissons$coeffs$D[1:4,1:4]))
hep <- mixed_model(fixed = hepatomegaly ~ drug * time,   # hepatomegaly
                   random = ~1|id,
                   data = PBC, family = binomial())
M <- list(lsb, alb, ast, plt, hep)                       # List of all responses

sdt <- PBC[!duplicated(PBC$id), ]                        # Survival data
ph <- coxph(Surv(survtime, status) ~ drug, sdt)          # PH fit
startup <- round(proc.time()[3] - pt, 3)                 # Elapsed time for startup.
# 205.897s

# Save everything before JMb2 fit
all.fits <- list(Gaussians = Gaussians,
                 Poissons = Poissons,
                 Binomials = Binomials,
                 Reduced = reduced.model,
                 JMb2 = list(
                   constituents = M,
                   ph = ph
                 ))
save(all.fits, file = '~/Downloads/prefit.RData')

jmb <- jm(ph, M, 'time', id_var = 'id',
          control = list(n_chains = 3L,
                         n_iter = 7000L,
                         n_burnin = 500L))

# SATURATED ---------------------------------------------------------------
saturated7 <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id),
       AST ~ drug * time  + (1 + time|id),
       prothrombin ~ drug * time  + (1 + time|id),
       platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1, 4))|id),
       alkaline ~ drug * time  + (1 + time|id),
       hepatomegaly ~ drug * time + (1|id)),
  surv.formula,
  PBC,
  family = list("gaussian", "gaussian", "gaussian", "gaussian",
                "poisson", "poisson", "binomial")
)


saturated8 <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id),
       AST ~ drug * time  + (1 + time|id),
       prothrombin ~ drug * time  + (1 + time|id),
       platelets ~ drug * ns(time, knots = c(1,4)) + (1 + ns(time, knots = c(1, 4))|id),
       alkaline ~ drug * time  + (1 + time|id),
       hepatomegaly ~ drug * time + (1|id),
       spiders ~ drug * time  + (1|id)),
  surv.formula,
  PBC,
  family = list("gaussian", "gaussian", "gaussian", "gaussian",
                "poisson", "poisson", "binomial", "binomial")
)