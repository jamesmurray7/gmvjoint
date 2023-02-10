# Investigating MVN-ness!
D <- diag(c(0.16, 0.04, 0.25, 0.05, 0.20, 0.04))
D[lower.tri(D, )] <- c(0.03, 0.02, 0.04, 0.0, 0.0, # Gaussian (Intercept) and others.
                       0.03, 0.00, -0.06,  0.0  ,  # Gaussian (slope) and others.
                       0.08, 0.05,   0.0,          # Poisson (Intercept) and others.
                       0.04, 0.0,0.0)                   # Poisson slope and binomial intercept

# Fixed effects
beta <- rbind(
  c(0.2, -0.2, 0.2, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(2, -0.1, 0.1, -0.2)
  # c(1, -1, 1, -1)
)

# Dispersion (only needed for sigma)
sigma <- c(0.16, 0, 0)

gamma <- c(-0.5, 0.25, 0.40)
zeta <- c(-0.0, 0.3)

family <- list('gaussian', 'poisson', 'poisson')

# Try for one...
load("~/Downloads/n_100_mi_6.RData")
bb <- survfit(Surv(dm$ph$survdata$survtime,dm$ph$survdata$status)~1)
get.post.b <- function(x){ # x a SINGLE joint object
  if(!inherits(x, 'joint')) stop("111")
  dm <- x$dmats
  sf <- survfit(Surv(dm$ph$survdata$survtime, dm$ph$survdata$status) ~ 1)
  #exp(-2.5 + 0.1 * dm$ph$ft))
  sv <- surv.mod(dm$ph, lapply(x$ModelInfo$long.formulas, parseFormula), diff(c(0,sf$cum[sf$n.ev>=1])))
  b <- lapply(1:x$ModelInfo$n, function(i) x$REs[i,])
  postb <- lapply(1:x$ModelInfo$n, function(i){
    optim(rep(0,sv$q), joint_density, joint_density_ddb, 
          Y = dm$long$Y[[i]], X = dm$long$X[[i]], Z = dm$long$Z[[i]],
          beta = c(utils::stack(as.data.frame(t(beta)))$values),
          D = D, sigma = list(0.16, 0, 0), family = family, 
          Delta = dm$ph$Delta[[i]], S = sv$S[[i]], Fi = sv$Fi[[i]], l0i = sv$l0i[[i]],
          SS = sv$SS[[i]], Fu = sv$Fu[[i]], haz = sv$l0u[[i]], gamma_rep = rep(gamma, c(2,2,2)),
          zeta = zeta[2], beta_inds = lapply(x$ModelInfo$inds$beta, function(x) x-1),
          b_inds = lapply(x$ModelInfo$inds$b, function(x) x-1),
          K = 3, method = 'BFGS', hessian = F)$par
  })
  out <- do.call(rbind, postb)
  colnames(out) <- paste0('b[',paste0(rep(1:3, each=2), rep(0:1, 2)),']')
  out
}

bs <- lapply(fits_n100m6, get.post.b)

library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_csda())

allbs <- as.data.frame(do.call(rbind, bs))
allbs$id <- rep(1:100, each = 100)

theor <- do.call(cbind,lapply(1:6, function(x){
  b <- allbs[,x,drop=T]
  den <- density(b)
  dnorm(seq(min(den$x), max(den$x), length.out = 1e3),
        mean = 0, sd = sqrt(D[x,x]))
})) %>% as.data.frame
names(theor) <- names(allbs)[1:6]
theorlong <- theor %>% pivot_longer(everything())

allbs %>% pivot_longer(-id) %>% arrange(id, name) %>% 
  filter(id <= 25) %>% 
  ggplot(aes(x = value, group = id)) + 
  geom_density(colour = alpha('black', 0.25)) + 
  geom_density(data = theorlong, mapping = aes(x = value, group = NULL),
               colour = 'red') + 
  facet_wrap(~name, scales = 'free', labeller = label_parsed, ncol = 2)










