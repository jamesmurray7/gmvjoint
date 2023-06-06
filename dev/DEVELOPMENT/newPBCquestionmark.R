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

# Relative conv crit to more params.
reduced.model.lower.tol <- joint(reduced.long,
                                 surv.formula, PBC, list("gaussian", "gaussian", "gaussian",
                                                         "poisson", "binomial"),
                                 control = list(tol.thr = 0.01))

summary(reduced.model)
xtable(reduced.model, vcov = TRUE, size = 'footnotesize',
       booktabs = FALSE)


# JMbayes2 fit ------------------------------------------------------------
library(JMbayes2)
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

jmb <- jm(ph, M, time_var = 'time', id_var = 'id',
          data_Surv = sdt,
          control = list(n_chains = 1L,
                         n_iter = 11000L,
                         n_burnin = 1000L))

save(jmb, file = '~/Downloads/New_jmb.RData')
load('~/Downloads/New_jmb.RData')

# JMbayes2:::summary.jm provides mean (SD) [95% CrI] for {beta, sigma, gamma, zeta}
# But we need to work out these for vech(D) separately.
sjmb <- summary(jmb)

# Start with D
jmb.D.mean <- jmb$statistics$Mean$D
jmb.D.SD <- jmb$statistics$SD$D
jmb.D.lb <- jmb$statistics$CI_low$D
jmb.D.ub <- jmb$statistics$CI_upp$D

# Get something xtable-ready.
.to3dp <- function(x) format(round(x, 3), nsmall = 3)
xt <- xtable(reduced.model, vcov = T)
chunks <- lapply(1:5, function(i){
  x <- reduced.model$ModelInfo$inds$R$b[[i]]
  # Start with D
  Dx <- sjmb$D[x,x]; vDx <- vech(Dx)
  inds <- which(vech(sjmb$D) %in% vDx)
  xSD <- jmb.D.SD[inds]
  xLB <- jmb.D.lb[inds]
  xUB <- jmb.D.ub[inds]
  Parameterjmb <- paste0("D_{", i,",",apply(which(lower.tri(Dx, T), arr.ind = T) - 1, 1, paste, collapse=''),"}")
  MSD <- paste0(.to3dp(vDx), " (", .to3dp(xSD), ")")
  CrI <- paste0('[', .to3dp(xLB),', ',.to3dp(xUB),']')
  
  JMbD <- cbind(Parameterjmb, MSD, CrI)
  # Now beta/sigma(^2)/gamma (in that order)
  gammas <- sjmb$Survival[i+1, c("Mean", "StDev", "2.5%", "97.5%"),drop=F]
  lookup <- paste0("Outcome",i)
  Out <- sjmb[[lookup]][,c("Mean", "StDev", "2.5%", "97.5%")]
  Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")] <- Out[grepl("sigma", row.names(Out)), c("Mean", "2.5%", "97.5%")]^2
  this <- rbind(Out, gammas)
  MSD <- paste0(.to3dp(this$Mean), ' (', .to3dp(this$StDev), ')')
  CrI <- paste0('[',.to3dp(this$`2.5%`), ', ', .to3dp(this$`97.5%`), ']')
  Parameterjmb <- row.names(this)
  JMbrest <- cbind(Parameterjmb, MSD, CrI)
  out <- cbind(xt$RespChunks[[i]], rbind(JMbD, JMbrest))
  if(i == 5){
    zeta <- sjmb$Survival[1,c("Mean", "StDev", "2.5%", "97.5%"),drop=F]
    MSD <- paste0(.to3dp(zeta$Mean), ' (', .to3dp(zeta$StDev), ')')
    CrI <- paste0('[',.to3dp(zeta$`2.5%`), ', ', .to3dp(zeta$`97.5%`), ']')
    out <- rbind(out, cbind(xt$zeta, Parameterjmb=row.names(zeta), MSD, CrI))
  }
  out
})

to.xt2 <- do.call(rbind, chunks)
to.xt2 <- to.xt2[,-4]
to.xt2[,1] <- paste0("$\\", to.xt2[,1], "$")
xt.to.xt <- xtable(to.xt2)
print(xt.to.xt, 
      include.rownames = FALSE,
      sanitize.text.function = identity)

# Summary param-by-param
srm <- summary(reduced.model)
for(i in 1:6){
  if(i == 6){
    cat("\nPress Enter for survival sub-model\n")
    xx <- readline()
    cat("Approximate EM:\n")
    print(srm$Survs)
    cat("\n\nJMbayes2:\n")
    print(sjmb$Survival)
    cat("\n")
  }else{
    cat(sprintf("\nPress enter for response %d", i))
    xx <- readline()
    cat("Approximate EM:\n")
    print(srm$Longits[[i]])
    cat("\n\nJMbayes2:\n")
    print(sjmb[[paste0("Outcome", i)]])
    cat("\n")
  }
}

library(ggplot2)
library(dplyr)




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




















# INLA -- DOESNT WORK -----------------------------------------------------
# library(INLAjoint)
# data(pbc2, package = 'JM')
# pbc2temp <- subset(pbc2, select = c("id", "year", "drug", "serBilir", 
#                                     "albumin", "SGOT", "platelets", "hepatomegaly",
#                                     "years", "status"))
# pbc2temp$serBilir <- log(pbc2temp$serBilir)
# pbc2temp$SGOT <- log(pbc2temp$SGOT)
# pbc2temp$status2 <- ifelse(pbc2temp$status == 'dead', 1, 0)
# pbc2temp <- na.omit(pbc2temp)
# # Survival data
# SD <- pbc2temp[!duplicated(pbc2temp$id), c("id", "years", "status2", "drug")]
# # time functoins
# f1 <- function(x) x^2
# NS <- splines::ns(pbc2temp$year, knots = c(1, 4))
# f2 <- function(x) predict(NS, x)[,1]
# f3 <- function(x) predict(NS, x)[,2]
# f4 <- function(x) predict(NS, x)[,3]
# 
# death <- inla.surv(time = c(SD$years), event = c(SD$status2))
# # Can't get this to work with SRE, and INLA crashes with CV (as below)...
# inj <- joint(
#   formSurv = death ~ drug,
#   formLong = list(
#     serBilir ~ drug * (1 + year + f1(year)) + (1 + year + f1(year) | id),
#     albumin ~ drug * (1 + year) + (1 + year | id),
#     SGOT ~ drug * (1 + year) + (1 + year | id),
#     platelets ~ drug * (1 + f2(year) + f3(year) + f4(year)) +
#       (1 + f2(year) + f3(year) + f4(year) | id),
#     hepatomegaly ~ drug * (1 + year) + (1 | id)
#   ), dataLong = pbc2temp, #dataSurv = SD,
#   id = 'id', timeVar = 'year', corLong = TRUE,
#   family = c("gaussian", "gaussian", "gaussian", "poisson", "binomial"),
#   link = rep("default", 5),
#   assoc = as.list(rep("CV", 5)), basRisk = "rw1",
#   control = list(int.strategy = 'eb', cfg = TRUE), NbasRisk = 15,
# )