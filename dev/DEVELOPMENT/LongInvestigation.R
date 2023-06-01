rm(list=ls())
data("PBC", package = 'gmvjoint')
library(glmmTMB)
PBC$AST <- log(PBC$SGOT)
# PBC$prothrombin <- (PBC$prothrombin * .1)^(-4)


# Binary ------------------------------------------------------------------
# No matter what, we fit binary with intercept only RE.

# Counts ------------------------------------------------------------------
# Platelets
P1 <- glmmTMB(platelets ~ drug * time + (1 + time|id), PBC, poisson)
VarCorr(P1); logLik(P1); AIC(P1); BIC(P1)
P2 <- glmmTMB(platelets ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, poisson)
VarCorr(P2); logLik(P2); AIC(P2); BIC(P2)
P3 <- glmmTMB(platelets ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, poisson)
VarCorr(P3); logLik(P3); AIC(P3); BIC(P3)
# Spline?
# Residuals always look good
par(mfrow=c(3,1)); plot(residuals(P1, 'pearson'));plot(residuals(P2, 'pearson'));plot(residuals(P3,'pearson'));par(mfrow=c(1,1))

# Alkaline
A1 <- glmmTMB(alkaline ~ drug * time + (1 + time|id), PBC, glmmTMB::nbinom2())
VarCorr(A1); logLik(A1); AIC(A1); BIC(A1)
A2 <- glmmTMB(alkaline ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, glmmTMB::nbinom2())
VarCorr(A2); logLik(A2); AIC(A2); BIC(A2)
A3 <- glmmTMB(alkaline ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, glmmTMB::nbinom2()) # Fails for Poisson
VarCorr(A3); logLik(A3); AIC(A3); BIC(A3)
par(mfrow=c(2,1)); plot(residuals(A1, 'pearson'));plot(residuals(A2, 'pearson'));par(mfrow=c(1,1))
# Negbin fits a lot better than Poisson; residuals here point towards dropping this response altogether...

# Continuous --------------------------------------------------------------
# serBilir
S1 <- glmmTMB(log(serBilir) ~ drug * time + (1 + time|id), PBC, gaussian)
VarCorr(S1); logLik(S1); AIC(S1); BIC(S1)
S2 <- glmmTMB(log(serBilir) ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, gaussian)
VarCorr(S2); logLik(S2); AIC(S2); BIC(S2)
S3 <- glmmTMB(log(serBilir) ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, gaussian)
VarCorr(S3); logLik(S3); AIC(S3); BIC(S3)
par(mfrow=c(3,1)); plot(residuals(S1, 'pearson'));plot(residuals(S2, 'pearson'));plot(residuals(S3,'pearson'));par(mfrow=c(1,1))
# Quadratic

# Aspartate aminotransfarase (log(AST))
SG1 <- glmmTMB(AST ~ drug * time + (1 + time|id), PBC, gaussian)
VarCorr(SG1); logLik(SG1); AIC(SG1); BIC(SG1)
SG2 <- glmmTMB(AST ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, gaussian)
VarCorr(SG2); logLik(SG2); AIC(SG2); BIC(SG2)
SG3 <- glmmTMB(AST ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, gaussian)
VarCorr(SG3); logLik(SG3); AIC(SG3); BIC(SG3)
par(mfrow=c(3,1)); plot(residuals(SG1, 'pearson'));plot(residuals(SG2, 'pearson'));plot(residuals(SG3,'pearson'));par(mfrow=c(1,1))
# Linear or quadratic

# Albumin
A1 <- glmmTMB(albumin ~ drug * time + (1 + time|id), PBC, gaussian)
VarCorr(A1); logLik(A1); AIC(A1); BIC(A1)
A2 <- glmmTMB(albumin ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, gaussian)
VarCorr(A2); logLik(A2); AIC(A2); BIC(A2)
A3 <- glmmTMB(albumin ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, gaussian) # Fails for Poisson
VarCorr(A3); logLik(A3); AIC(A3); BIC(A3)
par(mfrow=c(3,1)); plot(residuals(A1, 'pearson'));plot(residuals(A2, 'pearson'));plot(residuals(A3, 'pearson'));par(mfrow=c(1,1))
# Linear (rest fail!)

# Prothrombin
# No transform
P1 <- glmmTMB(prothrombin ~ drug * time + (1 + time|id), PBC, gaussian)
VarCorr(P1); logLik(P1); AIC(P1); BIC(P1)
P2 <- glmmTMB(prothrombin ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id), PBC, gaussian)
VarCorr(P2); logLik(P2); AIC(P2); BIC(P2)
P3 <- glmmTMB(prothrombin ~ drug * ns(time, 3) + (1 + ns(time, 3)|id), PBC, gaussian)
VarCorr(P3); logLik(P3); AIC(P3); BIC(P3)
# Residuals always look good
par(mfrow=c(3,1)); plot(residuals(P1, 'pearson'));plot(residuals(P2, 'pearson'));plot(residuals(P3,'pearson'));par(mfrow=c(1,1))


