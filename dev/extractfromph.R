beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
gamma <- c(0.3, -0.3)
D <- diag(c(0.25, 0.09, 0.25, 0.05))
family <- list('gaussian', 'poisson')
data <- simData(ntms = 10, beta = beta, D = D, n = 100,
                family = family, zeta = c(0, -0.2),
                sigma = c(0.16, 0), gamma = gamma)$data

# Specify formulae and target families
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id)  # Poisson
)

data$failed <- data$status
data$year <- data$survtime
data$survtime <- data$status <- NULL

surv.formula <- Surv(year, failed) ~ bin
# surv.formula <- Surv(survtime, status) ~ bin

fit <- joint(long.formulas, surv.formula, data, family)


# Now "Get back" to our survival input...
ph <- fit$dmats$ph$ph
extract.surv.process <- function(ph){
  if(!inherits(ph, 'coxph')) stop("ph must be object of class 'coxph'.")
  call <- deparse(ph$formula)
  call <- gsub('\\(|\\)', '',
               regmatches(call, regexpr("\\(.*\\)", call)))
  splitcall <- trimws(el(strsplit(call, '\\,')))
  Time <- splitcall[1]; Status <- splitcall[2]
  return(list(Time = Time, Status = Status))
}
