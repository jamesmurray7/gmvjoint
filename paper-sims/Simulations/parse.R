# Restate true parameter values
D.true <- diag(c(0.25, 0.09, 0.50, 0.10, 2.00))
D.true[1,3] <- D.true[3,1] <- D.true[1,5] <- D.true[5,1] <- D.true[3,5] <- D.true[5,3] <- 0.25
beta.true <- c(2, -0.1, 0.1, -0.2, 2, -0.1, 0.1, -0.2, 1, -1, 1, -1)
sigma.true <- .16
gamma.true <- c(.5,-.5,.5)
zeta.true <- -.2
target <- c(vech(D.true), beta.true, sigma.true, gamma.true, zeta.true)
target.mat <- t(apply(t(matrix(target)),2,rep,100))
# Parsing fits.
qz <- qnorm(.975)
# Function to extract estimates
extract.estimates <- function(L) sapply(L, function(x){
  setNames(c(vech(x$coeffs$D), x$coeffs$beta, unlist(x$coeffs$sigma)[1], x$coeffs$gamma, x$coeffs$zeta),
           names(x$SE))
})
# Function to extract SEs
extract.SE <- function(L) sapply(L, function(x){
  x$SE
})

ests <- lapply(fits, extract.estimates)
SEs <- lapply(fits, extract.SE)

# Empirical Means, SDs and average estimated SE.
emp.mean <- lapply(ests, rowMeans)
emp.sd <- lapply(ests, function(x) apply(x, 1, sd))
avg.se <- lapply(SEs, rowMeans)

# 95% coverage probability
mapply(function(estimate, se){
  lb <- estimate - qz * se; ub <- estimate + qz * se
  rowSums(lb <= target.mat & ub >= target.mat)/100
}, estimate = ests, se = SEs)
