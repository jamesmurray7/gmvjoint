load('/data/c0061461/ADNI.RData')

prepdata <- function(vars){
  this <- na.omit(ADNI[, c('id', 'survtime', 'month', 'time', 'status', 'gender', 'bin', vars, 'RID')])
  cn <- colnames(this)
  colnames(this)[length(cn)] <- 'originalRID'
  uids <- unique(this$id)
  key <- data.frame(id = uids, newid = 1:length(uids))
  this <- merge(this, key, 'id')
  this$id <- this$newid
  this$newid <- NULL
  this
}


# ADAS13 ------------------------------------------------------------------
set1 <- prepdata(c('ADAS13', 'PTEDUCAT', 'AGE'))
surv.formula <- Surv(survtime, status) ~ bin + gender + PTEDUCAT + AGE

long.formulaY <- list(
  Y ~ time + bin + gender + PTEDUCAT + AGE + (1 + time|id)
)
set1$Y <- c(scale(set1$ADAS13))

joint.G <- joint(long.formulaY, surv.formula, set1, list("gaussian"), control = list(verbose = T))

# Now treat it as a count.
long.formula <- list(
  ADAS13 ~ time + bin + gender + PTEDUCAT + AGE + (1 + time|id)
)

joint.P <- joint(long.formula, surv.formula, set1, list("poisson"), control = list(verbose = T))

joint.GP <- joint(long.formula, surv.formula, set1, list("genpois"), control = list(verbose = T))

a <- bootAUCdiff(list(joint.G, joint.GP), set1, 3, 3, nboot = 100)
b <- bootAUCdiff(list(joint.G, joint.GP), set1, 5, 2, nboot = 100)
d <- bootAUCdiff(list(joint.G, joint.GP), set1, 6, 2, nboot = 100)

# Any different with less impactful predictors removed from model?
joint.G2 <- joint(list(
  Y ~ time + bin + PTEDUCAT + (1 + time|id)
), Surv(survtime, status) ~ bin + AGE, data = set1, family = list("gaussian"))

joint.P2 <- joint(list(
  ADAS13 ~ time + bin + PTEDUCAT + (1 + time|id)
), Surv(survtime, status) ~ bin + AGE, data = set1, list("poisson"), control = list(verbose = F))

joint.GP2 <- joint(list( # Var[Y] is approximately 70% of E[Y].
  ADAS13 ~ time + bin + PTEDUCAT + (1 + time|id)
), Surv(survtime, status) ~ bin + AGE, data = set1, list("genpois"), control = list(verbose = T))


# MMSE --------------------------------------------------------------------
# bounded [0,30]
set2 <- prepdata(c("MMSE", "PTEDUCAT", "AGE"))

set2$Y <- c(scale(set2$MMSE))
set2$Lost <- c(max(set2$MMSE) - set2$MMSE)

long.formula.Y.saturated <- list(
  Y ~ time + bin + gender + PTEDUCAT + AGE + (1 + time|id)
)
long.formula.M.saturated <- list(
  Lost ~ time + bin + gender + PTEDUCAT + AGE + (1+time|id)
)

surv.saturated <- Surv(survtime, status) ~ bin + gender + PTEDUCAT + AGE

joint.G <-  joint(long.formula.Y.saturated, surv.saturated, data = set2, family = list('gaussian'))
(joint.P <-  joint(long.formula.M.saturated, surv.saturated, data = set2, family = list('poisson')))
(joint.GP <- joint(long.formula.M.saturated, surv.saturated, data = set2, family = list('genpois')))

a <- bootAUCdiff(fits = list(joint.G, joint.GP), data = set2, Tstart = 3, delta = 3, nboot = 25)
# GP better here, but we need to remodel the response which probably isn't good!



# RAVLT -------------------------------------------------------------------
set3 <- prepdata(c("PTEDUCAT", "AGE", "RAVLT.immediate",
                   "RAVLT.learning", "RAVLT.forgetting", "RAVLT.perc.forgetting"))
set3

joint.P <- joint(
  list(RAVLT.immediate ~ time + PTEDUCAT + AGE + bin + (1 + time|id)),
  Surv(survtime, status) ~ bin + AGE,
  data = set3, family = list("poisson")
)

joint.GP <- joint(
  list(RAVLT.immediate ~ time + PTEDUCAT + AGE + bin + (1 + time|id)),
  Surv(survtime, status) ~ bin + AGE,
  data = set3, family = list("genpois")
)

set3$Y.imm <- c(scale(set3$RAVLT.immediate))
set3$Y.perc
joint.G <- joint(
  list(Y.imm ~ time + PTEDUCAT + AGE + bin + (1 + time|id)),
  Surv(survtime, status) ~ bin + AGE,
  data = set3, family = list("Gamma")
)

a <- bootAUCdiff(list(joint.P, joint.GP), set3, 6, 2) # Uniformally better to model with GP.






