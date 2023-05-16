rm(list=ls())
library(splines)
data(PBC)

PBC <- na.omit(PBC[,c("id", "survtime", "status", "drug", "sex", "age", "time", 'ascites',
                      "hepatomegaly", "spiders", "serBilir",
                      "albumin", "alkaline", "SGOT", "platelets", "prothrombin")])
PBC$serBilir <- log(PBC$serBilir)
PBC$prothrombin <- (PBC$prothrombin * .1)^ (-4)
PBC$AST <- log(PBC$SGOT)

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

PBC$ascites <- NULL
PBC <- na.omit(PBC)

Gaussian.long.formulas <- list(
  serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
  albumin ~ drug * time + (1 + time|id),
  prothrombin ~ drug * ns(time, 3) + (1 + ns(time, 3)|id),
  AST ~ drug * time  + (1 + time|id)
)

Poisson.long.formulas <- list(
  platelets ~ drug * time  + (1 + time|id),
  alkaline ~ drug * time  + (1 + time|id)
)

Binomial.long.formulas <- list(
  hepatomegaly ~ drug * time  + (1|id),
  spiders ~ drug * time  + (1|id)
)

surv.formula <- Surv(survtime, status) ~ drug
control <- list(verbose=T)

Gaussians <- joint(Gaussian.long.formulas, surv.formula, PBC, 
                   list("gaussian", "gaussian", "gaussian", "gaussian"))

Poissons <- joint(Poisson.long.formulas,
                  surv.formula, PBC, list("poisson", "poisson"))

Binomials <- joint(Binomial.long.formulas, surv.formula,
                   PBC, list("binomial", "binomial"))

source('testing/fit2xtab.R')

fit2xtab(Gaussians, 16)
fit2xtab(Poissons)
fit2xtab(Binomials)

# Take Gaussian: {serBilir, albumin}, Poisson: {platelets} and Binomial {hepatomegaly} forward

reduced.model <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id),
       platelets ~ drug * time  + (1 + time|id),
       hepatomegaly ~ drug * time  + (1|id)),
  surv.formula, PBC, list("gaussian", "gaussian", "poisson", "binomial"))

summary(reduced.model)
fit2xtab(reduced.model, )
test <- cond.ranefs(reduced.model,tune=0.80)

final.biv.model <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id)),
  surv.formula, PBC, list("gaussian", "gaussian"), control = control
)
summary(final.biv.model)
test <- cond.ranefs(final.biv.model,tune=0.80)


# Full seven-variate? -----------------------------------------------------

all.long.formulas <- c(Gaussian.long.formulas, Poisson.long.formulas, Binomial.long.formulas[[1]])

all.fit <- joint(all.long.formulas,
                 surv.formula, PBC, list("gaussian", "gaussian", "gaussian", "gaussian",
                                         "poisson", "poisson", "binomial"))
save(all.fit, file = '/data/c0061461/GLMM_Paper_Sims/Revision2/PBCallfits.RData')