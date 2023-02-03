all.sims <- expand.grid(mi = c(5,10,15), 
                        n =  c(100, 250, 500))
            # G ------> P ------->  B

# Families
family <- list("gaussian", "poisson", "binomial")

# VarCorr matrix for REs
D <- diag(c(0.25, 0.04, 0.25, 0.04, 0.25))
D[lower.tri(D, )] <- c(0.05, 0.04, -0.01, 0.0,  # Gaussian (Intercept) and others.
                       0.00, -0.03, 0.0,        # Gaussian (slope) and others.
                       -0.02, 0.1,              # Poisson (Intercept) and others.
                       0)                       # Poisson slope and binomial intercept

D[upper.tri(D)] <- t(D)[upper.tri(D)]


D <- diag(5)
D[1, 1] <- D[3, 3] <- D[5, 5] <- 0.5^2
D[2, 2] <- D[4, 4] <- 0.2^2
D[1, 3] <- D[3, 1] <- -0.5 * 0.5 * 0.5
D[1, 5] <- D[5, 1] <- 0.5^3
D[3, 5] <- D[5, 3] <- -0.5*(0.5^2)



# Check this ok
isSymmetric(D) && all(eigen(D)$val>0) && det(D) > 0

# Fixed effects
beta <- rbind(
  c(0.2, -0.2, 0.2, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(1, -1, 1, -1)
)


beta <- rbind(c(0, 1, 1, 1),
              c(0, -1, 0, 0.5),
              c(0, 0.1, 0.5, -0.5))

# Dispersion (only needed for sigma)
sigma <- c(0.16, 0, 0)

# gamma & zeta
gamma <- c(-0.5, 0.75, -0.50)
zeta <- c(-0.0, 0.3)

# Specify random effects formulae (W_k)
random.formula <- list(~time, ~time, ~1)

test <- simData(family = family, sigma = sigma, ntms = 15, n = 100,
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
a <- fitfn(test$data,control = list(correlated=F))

















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




