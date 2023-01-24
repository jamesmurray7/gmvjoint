library(gmvjoint)
data(PBC)
PBC <- na.omit(PBC[,c("id", "survtime", "drug", "sex", "time", "ascites",
                      "spiders", "serBilir", "albumin", "alkaline", "SGOT",
                      "platelets", "prothrombin", "hepatomegaly", "status", "age")])
PBC$id <- as.numeric(as.character(PBC$id))

PBC$serBilir <- log(PBC$serBilir)
PBC$SGOT <- log(PBC$SGOT)
PBC$prothrombin <- (0.1 * PBC$prothrombin)^(-4)

library(splines)

make.formula <- function(response, type){
  f <- switch(type,
              linear = paste0(response, '~ drug * time + (1 + time|id)'),
              spline = paste0(response, '~ drug * ns(time, 3) + (1 + ns(time, 3)|id)'),
              quad   = paste0(response, '~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id)'),
              int    = paste0(response, '~ drug * time + (1|id)'))
  as.formula(f)  
}

# Longits -----------------------------------------------------------------
# Conts
SB <- glmmTMB(make.formula('serBilir', 'quad'), PBC)
SG <- glmmTMB(make.formula('SGOT', 'linear'), PBC)
alb <- glmmTMB(make.formula('albumin', 'linear'), PBC)
pro <- glmmTMB(make.formula('prothrombin', 'spline'), PBC)
# Counts
pla <- glmmTMB(make.formula('platelets', 'linear'), PBC, family = poisson)
alk <- glmmTMB(make.formula('alkaline', 'linear'), PBC, family = poisson)
# Bins
spi <- glmmTMB(make.formula('spiders', 'int'), PBC, family = binomial)
asc <- glmmTMB(make.formula('ascites', 'int'), PBC, family = binomial)
hep <- glmmTMB(make.formula('hepatomegaly', 'int'), PBC, family = binomial)

extract.ranef <- function(x){
  out <- as.matrix(ranef(x)$cond$id)
  colnames(out) <- paste0(names(x$modelInfo$respCol), '_', 
                          colnames(out))
  out
}

conts <- list(SB = extract.ranef(SB), SG = extract.ranef(SG), alb = extract.ranef(alb), 
              pro = extract.ranef(pro))
conts.corr <- round(cov2cor(cov(do.call(cbind, conts))), 4)   # |Largest| appx. .4

counts <- list(pla = extract.ranef(pla), alk = extract.ranef(alk)) 
counts.corr <- round(cov2cor(cov(do.call(cbind, counts))), 4) # Small

bins <- list(spi = extract.ranef(spi), asc = extract.ranef(asc),
             hep = extract.ranef(hep)) 
bins.corr <- round(cov2cor(cov(do.call(cbind, bins))), 4) # Moderate


# JMs ---------------------------------------------------------------------
# bins all modelled together
all.bin <- joint(
  long.formulas = list(
    make.formula('spiders', 'int'),
    make.formula('ascites', 'int'),
    make.formula('hepatomegaly', 'int')
  ),
  surv.formula = Surv(survtime, status) ~ drug,
  data = PBC, #control = list(verbose = T),
  family = list('binomial', 'binomial', 'binomial')
)
summary(all.bin) # Take Ascites.
cov2cor(all.bin$coeffs$D)


# conts all modelled together.
all.cont <- joint(
  long.formulas = list(
    make.formula('serBilir', 'quad'),
    make.formula('SGOT', 'linear'),
    make.formula('albumin', 'linear'),
    make.formula('prothrombin', 'spline')
  ),
  surv.formula = Surv(survtime, status) ~ drug,
  data = PBC, #control = list(verbose = T),
  family = list('gaussian', 'gaussian', 'gaussian', 'gaussian')
)
summary(all.cont)
cov2cor(all.cont$coeffs$D)
# Would keep all; remove SGOT.

# counts all modelled
all.count <- joint(
  long.formulas = list(
    make.formula('platelets', 'linear'),
    make.formula('alkaline', 'linear')
  ),
  surv.formula = Surv(survtime, status) ~ drug,
  data = PBC, control = list(verbose = T),
  family = list('poisson', 'poisson')
)
summary(all.count) # Take platelets
cov2cor(all.count$coeffs$D)

# candidate 1
fit1 <- joint(
  long.formulas = list(
    make.formula('serBilir', 'quad'),
    make.formula('albumin', 'linear'),
    make.formula('prothrombin', 'spline'),
    make.formula('platelets', 'linear'),
    make.formula('ascites', 'int')
  ), surv.formula = Surv(survtime, status) ~ drug,
  data = PBC, control = list(verbose = T),
  family = list('gaussian', 'gaussian', 'gaussian', 'poisson',
                'binomial')
)
summary(fit1)
# Take serum bilirubin, albumin and ascites

# candidate 2 // reduced model
# An "update.joint" function could be nice?
fit2 <-  joint(
  long.formulas = list(
    make.formula('serBilir', 'quad'),
    make.formula('albumin', 'linear'),
    make.formula('ascites', 'int')
  ), surv.formula = Surv(survtime, status) ~ drug,
  data = PBC, control = list(verbose = F),
  family = list('gaussian', 'gaussian', 'binomial')
)
summary(fit2) # which is the _final_ model.



# Univariates only
fit <- joint(list(make.formula('SGOT', 'linear')), Surv(survtime, status)~drug,
      PBC,list('gaussian'))
summary(fit)
