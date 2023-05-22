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

# 4 variate Gaussian
Gaussians <- joint(Gaussian.long.formulas, surv.formula, PBC, 
                   list("gaussian", "gaussian", "gaussian", "gaussian"))

xtable(Gaussians, vcov = TRUE, max.row = 18, size = "tiny", booktabs = FALSE)
# Take forward: 

Poissons <- joint(Poisson.long.formulas,
                  surv.formula, PBC, list("poisson", "poisson"))

xtable(Poissons, vcov = TRUE, size = 'footnotesize', max.row = 9, booktabs = FALSE)

Binomials <- joint(Binomial.long.formulas, surv.formula,
                   PBC, list("binomial", "binomial"))
xtable(Binomials, vcov = TRUE, size = 'footnotesize', max.row = 7, booktabs = FALSE)

# Take Gaussian: {serBilir, albumin, prothrombin}, 
# Poisson: {platelets} and Binomial {hepatomegaly} forward...
reduced.long <- list(
     serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
     albumin ~ drug * time + (1 + time|id),
     prothrombin ~ drug * ns(time, 3) + (1 + ns(time, 3)|id),
     platelets ~ drug * time  + (1 + time|id),
     hepatomegaly ~ drug * time  + (1|id)
)

reduced.model <- joint(reduced.long,
  surv.formula, PBC, list("gaussian", "gaussian", "gaussian",
                          "poisson", "binomial"))

xtable(reduced.model, vcov = TRUE, size = 'footnotesize',
       booktabs = FALSE)

# JMbayes2 fit on final model ---------------------------------------------
library(JMbayes2)


pt <- proc.time()[3]                                     # Start timing
lsb <- lme(fixed = serBilir ~ drug * (time + I(time^2)), # log(serBilir)
           random = ~ time + I(time^2)|id,
           data = PBC)
alb <- lme(fixed = albumin ~ drug * time,                # albumin
           random = ~ time|id,
           data = PBC)
ptt <- lme(fixed = prothrombin ~ drug * (ns(time, 3)),   # prothrombin
           random = ~ 1 +  ns(time, 3) | id,
           data = PBC)
plt <- mixed_model(fixed = platelets ~ drug * time,      # platetlet count
                   random = ~time|id,                    # mixed_model struggles to fit for absolutely no reason, 
                   data = PBC, family = poisson(),       # so start at MLEs
                   initial_values = 
                     list(betas = c(5.52235114, -0.07439889, -0.06151984, -0.03921149),
                          D = matrix(c(0.146932, -0.00264, -0.00264, 0.022055), 2, 2)))
hep <- mixed_model(fixed = hepatomegaly ~ drug * time,   # hepatomegaly
                   random = ~1|id,
                   data = PBC, family = binomial())
M <- list(lsb, alb, ptt, plt, hep)                       # List of all responses

sdt <- PBC[!duplicated(PBC$id), ]                        # Survival data
ph <- coxph(Surv(survtime, status) ~ drug, sdt)          # PH fit
startup <- round(proc.time()[3] - pt, 3)                 # Elapsed time for startup.
# 57.023s (!!!!!!), 

# Fitting on default doesn't produce good mixing=
# jmb <- jm(ph, M, data_Surv = sdt, id_var = 'id',
#           time_var = 'time')
# Try for an arbitrarily long amount of iterations
jmb <- jm(ph, M, data_Surv = sdt, id_var = 'id',
          time_var = 'time',
          control = list(n_iter = 10000L,
                         n_burnin = 1000L))  
save(jmb, file = '/data/c0061461/GLMM_Paper_Sims/Revision2/ReducedModelJMb.RData')
# Mixing satisfactory for everything _but_ ns(time, 3)2 and ns(time, 3)3 on ptt.
# Already ran for 10,000 iterations, so just note this somewhere.
# Attach these to `reduced.model` above.
load("/data/c0061461/GLMM_Paper_Sims/Revision2/ReducedModelJMb.RData")
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
  x <- reduced.model$ModelInfo$ind$b[[i]]
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

to.xt <- do.call(rbind, chunks)
to.xt <- to.xt[,-4]
to.xt[,1] <- paste0("$\\", to.xt[,1], "$")
xt.to.xt <- xtable(to.xt)
print(xt.to.xt, 
      include.rownames = FALSE,
      sanitize.text.function = identity)

# Computation time
jmb$running_time[3] + startup
summary(reduced.model)$et[1] + summary(reduced.model)$et[2]
summary(reduced.model)$et[3]

# Full seven-variate? -----------------------------------------------------
# Keep spiders removed; it greatly increases computation time!!
all.long.formulas <- c(Gaussian.long.formulas, Poisson.long.formulas, Binomial.long.formulas)
all.long.formulas
# About 6 minutes
all.fit <- joint(all.long.formulas,
                 surv.formula, PBC, list("gaussian", "gaussian", "gaussian", "gaussian",
                                         "poisson", "poisson", "binomial"))
save(all.fit, file = '/data/c0061461/GLMM_Paper_Sims/Revision2/PBCallfits.RData')
xtable(all.fit, vcov = TRUE, booktabs = FALSE, max.row = 20)

