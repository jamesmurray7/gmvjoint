rm(list=ls())
library(splines)
library(xtable)
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
xtable(Gaussians)
xtable(Gaussians, vcov = T)

Poissons <- joint(Poisson.long.formulas,
                  surv.formula, PBC, list("poisson", "poisson"))
xtable(Poissons)

Binomials <- joint(Binomial.long.formulas, surv.formula,
                   PBC, list("binomial", "binomial"))
xtable(Binomials, max.row = 6)

# Take Gaussian: {serBilir, albumin}, Poisson: {platelets} and Binomial {hepatomegaly} forward

reduced.model <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id),
       platelets ~ drug * time  + (1 + time|id),
       hepatomegaly ~ drug * time  + (1|id)),
  surv.formula, PBC, list("gaussian", "gaussian", "poisson", "binomial"))

xtable(reduced.model)

final.biv.model <- joint(
  list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
       albumin ~ drug * time + (1 + time|id)),
  surv.formula, PBC, list("gaussian", "gaussian"), control = control
)
summary(final.biv.model)
xtable(final.biv.model)

# Full seven-variate? -----------------------------------------------------

all.long.formulas <- c(Gaussian.long.formulas, Poisson.long.formulas, Binomial.long.formulas)
# About 6 minutes
all.fit <- joint(all.long.formulas,
                 surv.formula, PBC, list("gaussian", "gaussian", "gaussian", "gaussian",
                                         "poisson", "poisson", "binomial", "binomial"))
save(all.fit, file = '/data/c0061461/GLMM_Paper_Sims/Revision2/PBCallfits.RData')
xtable(all.fit)
xtable(all.fit, max.row = 16)


# JMbayes2 final ----------------------------------------------------------
library(JMbayes2)

# log(serBilir)
lsb <- lme(fixed = serBilir ~ drug * (time + I(time^2)),
           random = ~ time + I(time^2)|id,
           data = PBC)
alb <- lme(fixed = albumin ~ drug * time,
           random = ~ time|id,
           data = PBC)
plt <- mixed_model(fixed = platelets ~ drug * time,
                   random = ~time|id,
                   data = PBC, family = poisson())
hep <- mixed_model(fixed = hepatomegaly ~ drug * time,
              random = ~1|id,
              data = PBC, family = binomial())
M <- list(lsb, alb, plt, hep)

sdt <- PBC[!duplicated(PBC$id), ]
ph <- coxph(Surv(survtime, status) ~ drug, sdt)

# Rhat < 1.05 for everything _but_ quadratic time term.
# (3.6 mins)
# jmb <- jm(ph, M, data_Surv = sdt, id_var = 'id',
#           time_var = 'time')
# Try for a bit longer?
jmb <- jm(ph, M, data_Surv = sdt, id_var = 'id',
          time_var = 'time',
          control = list(n_iter = 7000L))  #3.9 mins
(sjmb <- summary(jmb))

.ff <- function(x) format(round(x, 3), nsmall = 3)
fn <- function(x){
  x[grepl("sigma", rownames(x)), c("Mean", "2.5%", "97.5%")] <- x[grepl("sigma", rownames(x)), c("Mean", "2.5%", "97.5%")]^2
  `Mean (SD)` <- paste0(.ff(x$Mean), ' (', .ff(x$StDev), ')')
  CI <- paste0('[', .ff(x$`2.5%`), ', ', .ff(x$`97.5%`), ']')
  df <- data.frame(`Mean (SD)` = `Mean (SD)`, CI = CI,
                   row.names = rownames(x), stringsAsFactors = F)
  xt <- xtable(df)
  print(xt)
  xx <- readline('hit enter')
}
lapply(list(sjmb$Outcome1, sjmb$Outcome2, sjmb$Outcome3, sjmb$Outcome4), fn)
fn(sjmb$Survival)

