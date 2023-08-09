# Make a variable based on disease prog. at baseline ----------------------
data <- PBC
# Convert id to factor, or ordering gets altered => joint doesn't work properly.
data$id <- as.numeric(as.character(data$id)) 
hist.bl <- data[!duplicated(data$id), c("id", "histologic")]
names(hist.bl)[2] <- "BLdisease"
# Ensure this is a factor and merge back on
hist.bl$BLdisease <- factor(hist.bl$BLdisease, 1:4) # I assume this is ordinal 
data <- merge(data, hist.bl, "id")

long.formulas <- list(
  albumin ~ time * drug + edema + (1 + time|id) # factor in here to for parity's sake.
)
surv.formula <- Surv(survtime,status) ~ BLdisease + drug # factor + binary
control <- list(verbose=T)

m1 <- joint(long.formulas = long.formulas, surv.formula = surv.formula, 
            data = data, family = list("gaussian"), control = list(verbose = T))

# S3 methods appear to work decently well ->
summary(m1)
fixef(m1, "surv")

# Does this returns similar result to joineRML? ---------------------------
library(joineRML)

qq <- mjoint(
  formLongFixed = list(
    albumin ~ time * drug + edema
  ),
  formLongRandom = list(
    ~time|id
  ),
  formSurv = Surv(survtime,status) ~ BLdisease + drug,
  data = data,  timeVar = 'time',
  control = list(
    type = "sobol", tol2= 1e-2
  )
)

summary(qq)
summary(m1) # Looks like it
