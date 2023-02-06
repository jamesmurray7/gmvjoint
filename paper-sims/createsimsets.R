all.sims <- expand.grid(mi = c(5,10,15), 
                        n =  c(100, 250, 500))
            # G ------> P ------->  B

# Families
family <- list("gaussian", "poisson", "binomial")

# VarCorr matrix for REs
D <- diag(c(0.16, 0.04, 0.25, 0.05, 0.25))
D[lower.tri(D, )] <- c(0.03, 0.02, 0.04, 0.0,  # Gaussian (Intercept) and others.
                       0.03, 0.00, -0.06,        # Gaussian (slope) and others.
                       0.08, 0.05,              # Poisson (Intercept) and others.
                       0.04)                       # Poisson slope and binomial intercept

D[upper.tri(D)] <- t(D)[upper.tri(D)]

# Check this ok
isSymmetric(D) && all(eigen(D)$val>0) && det(D) > 0

# Fixed effects
beta <- rbind(
  c(0.2, -0.2, 0.2, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(1, -1, 1, -1)
)

# Dispersion (only needed for sigma)
sigma <- c(0.16, 0, 0)

# gamma & zeta
gamma <- c(-0.5, 0.25, 0.40)
zeta <- c(-0.0, 0.3)

# Specify random effects formulae (W_k)
random.formula <- list(~time, ~time, ~1)

test <- simData(family = family, sigma = sigma, ntms = 15, n = 250,
                beta = beta, D = D, gamma = gamma, zeta = zeta, random.formula = random.formula,
                theta = c(-2.5, 0.1))

# Function to fit ---------------------------------------------------------
fitfn <- function(d, ...){
  joint(list(Y.1 ~ time + cont + bin + (1 + time|id),
             Y.2 ~ time + cont + bin + (1 + time|id),
             Y.3 ~ time + cont + bin + (1|id)),
        Surv(survtime, status) ~ bin,
        d, family = list("gaussian", "poisson", "binomial"), ...)
}
a <- fitfn(test$data,control = list(verbose=T))


# Fit same using JMbayes2
data <- test$data
survdata <- test$surv.data
ph <- coxph(Surv(survtime, status) ~ bin, survdata)

gs <- lme(fixed=Y.1 ~ time + cont + bin,
          random = ~time|id, data = data, method = 'ML')
ps <- mixed_model(Y.2~time+cont+bin, data = data,
                  random = ~time|id, family = poisson())
bs <- mixed_model(Y.3~time+cont+bin, data = data,
                  random = ~1|id, family = binomial())


M <- list(gs, ps, bs)

jmb <- jm(ph, M, time_var = 'time')




# ignore, testing for MVN stuff. ------------------------------------------
dm <- a$dmats
bb <- survfit(Surv(dm$ph$survdata$survtime,dm$ph$survdata$status)~1)
sv <- surv.mod(dm$ph, lapply(a$ModelInfo$long.formulas, parseFormula), diff(c(0,bb$cum[bb$n.ev>=1])))#a$hazard[,2])#exp(-2.5 + 0.1 * dm$ph$ft))
b <- lapply(1:a$ModelInfo$n, function(i) a$REs[i,])
postb <- lapply(1:a$ModelInfo$n, function(i){
  print(i)
  print(joint_density(rep(0,sv$q), Y = dm$long$Y[[i]], X = dm$long$X[[i]], Z = dm$long$Z[[i]],
                      beta = c(utils::stack(as.data.frame(t(beta)))$values),
                      D = D, sigma = list(0.16, 0, 0), family = family, 
                      Delta = dm$ph$Delta[[i]], S = sv$S[[i]], Fi = sv$Fi[[i]], l0i = sv$l0i[[i]],
                      SS = sv$SS[[i]], Fu = sv$Fu[[i]], haz = sv$l0u[[i]], gamma_rep = rep(gamma, c(2,2,1)),
                      zeta = zeta[2], beta_inds = lapply(a$ModelInfo$inds$beta, function(x) x-1),
                      b_inds = lapply(a$ModelInfo$inds$b, function(x) x-1),
                      K = 3))
  optim(rep(0,sv$q), joint_density, joint_density_ddb, 
        Y = dm$long$Y[[i]], X = dm$long$X[[i]], Z = dm$long$Z[[i]],
        beta = c(utils::stack(as.data.frame(t(beta)))$values),
        D = D, sigma = list(0.16, 0, 0), family = family, 
        Delta = dm$ph$Delta[[i]], S = sv$S[[i]], Fi = sv$Fi[[i]], l0i = sv$l0i[[i]],
        SS = sv$SS[[i]], Fu = sv$Fu[[i]], haz = sv$l0u[[i]], gamma_rep = rep(gamma, c(2,2,1)),
        zeta = zeta[2], beta_inds = lapply(a$ModelInfo$inds$beta, function(x) x-1),
        b_inds = lapply(a$ModelInfo$inds$b, function(x) x-1),
        K = 3, method = 'BFGS')$par
})
allb <- do.call(rbind, postb)
par(mfrow=c(3,2))
for(q in 1:5){
     dens <- density(allb[,q])#, adjust = 1.5)
     plot(dens$x, dens$y, main = bquote("b"[.(q)]), ylab = 'density', 
                   xlab = 'posterior mode', type = 'l')
     curve(dnorm(x, mean =0, sd = sqrt(D[q,q])), add = T, col = 'red', lty = 3,
                     from = min(dens$x) , to = max(dens$x))
   }
par(mfrow=c(1,1))




