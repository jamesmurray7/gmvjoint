source('zzz/theme_csda.R')
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
theme_set(theme_csda())


# Model fit ---------------------------------------------------------------
ests <- lapply(fitsets, function(y) sapply(y, function(x){
  # s <- summary(x)
  co <- x$coeffs
  c(vech(co$D), c(co$beta), co$sigma[[1]], co$gamma, co$zeta)
}))

# SEs <- lapply(fitsets, function(y) sapply(y, function(x) sqrt(diag(solve(x$Itil)))))
SEs <- lapply(fitsets, function(y) sapply(y, function(x) x$SE))

# Lower and Upper bounds for CP
lb <- sapply(1:length(SEs), function(x){
  if(x %in% 4:6)
    ests[[x]] <<- ests[[x]][-8,]
  if(x %in% 7:9)
    ests[[x]] <<- ests[[x]][-6,]
  ests[[x]] - qnorm(.975) * SEs[[x]]
})

ub <- sapply(1:length(SEs), function(x){
  ests[[x]] + qnorm(.975) * SEs[[x]]
})

fams <- rep(c('G', "P", "B"), each = 3)

CP <- sapply(1:length(lb), function(x){
  fx <- fams[x]
  if(fx == 'B') vD <- 0.25 else vD <- c(0.25, 0, 0.05)
  if(fx != 'P') gamma <- -0.5 else gamma <- 0.5
  if(fx == 'G') sigma <- 0.16 else sigma <- NULL
  target <- c(vD, c(2, -0.1, 0.1, 0.2), sigma, gamma, 0.3)
  tmat <- matrix(target, nc = 100, nr = length(target), byr=F)
  lb[[x]] <= tmat & ub[[x]] >= tmat
})

sapply(CP, rowSums)


# "True" Conditional distns -----------------------------------------------

# Try for one...
get.post.b <- function(x){ # x a SINGLE joint object
  if(!inherits(x, 'joint')) stop("111")
  ff <- x$ModelInfo$family[[1]]
  if(ff != 'poisson') true.gamma <- -0.5 else true.gamma <- 0.5
  if(ff != 'binomial'){
    gamma.rep <- rep(true.gamma, 2)
    D <- diag(c(.25, .05))
  }else{
    gamma.rep <- true.gamma
    D <- matrix(.25, 1, 1)
  }
  dm <- x$dmats
  l0 <- exp(-2.5 + 0.1 * dm$ph$ft)
  sv <- surv.mod(dm$ph, lapply(x$ModelInfo$long.formulas, parseFormula), l0)
  b <- lapply(1:x$ModelInfo$n, function(i) x$REs[i,,drop=F])
  postb <- lapply(1:x$ModelInfo$n, function(i){
    optim(b[[i]], joint_density, joint_density_ddb, 
          Y = dm$long$Y[[i]], X = dm$long$X[[i]], Z = dm$long$Z[[i]],
          beta = c(2, -0.1, 0.1, 0.2),
          D = D, sigma = list(0.16), family = as.list(ff), 
          Delta = dm$ph$Delta[[i]], S = sv$S[[i]], Fi = sv$Fi[[i]], l0i = sv$l0i[[i]],
          SS = sv$SS[[i]], Fu = sv$Fu[[i]], haz = sv$l0u[[i]], gamma_rep = gamma.rep,
          zeta = 0.3, beta_inds = lapply(x$ModelInfo$inds$beta, function(x) x-1),
          b_inds = lapply(x$ModelInfo$inds$b, function(x) x-1),
          K = 1, method = 'BFGS', hessian = F,
          control = list(reltol = 1e-4))$par
  })
  out <- do.call(rbind, postb)
  if(ff != 'binomial') 
    colnames(out) <- c('b[0]', 'b[1]')
  else
    colnames(out) <- c('b[0]')
  out
}

# Gaussians ---------------------------------------------------------------
Dg <- diag(c(.25, .05))
bs5 <- lapply(fitsets$`n = 250, mi = 5, fam = gaussian`, get.post.b)
bs10 <- lapply(fitsets$`n = 250, mi = 10, fam = gaussian`, get.post.b)
bs15 <- lapply(fitsets$`n = 250, mi = 15, fam = gaussian`, get.post.b)

allbs5 <- as.data.frame(do.call(rbind, bs5))
allbs10 <- as.data.frame(do.call(rbind, bs10))
allbs15 <- as.data.frame(do.call(rbind, bs15))
allbs <- rbind(allbs5, allbs10, allbs15)
allbs$id <- rep(1:300, each = 250)
allbs$mi <- rep(c(5, 10, 15), each = (250 * 100))

theor <- do.call(rbind,sapply(1:2, function(x){
  b <- allbs[,x,drop=T]
  den <- density(b)
  value <- seq(min(den$x), max(den$x), length.out = 1e3)
  theory <- dnorm(value, mean = 0, sd = sqrt(Dg[x,x]))
  out <- data.frame(value = value, theory = theory)
  names(out)[1] <- colnames(allbs)[x]
  pivot_longer(out, -theory)
}, simplify = F))

allbs %>% 
  pivot_longer(`b[0]`:`b[1]`) %>%
  # filter(id %in% c(1:10, 100:110, 200:210)) %>%
  ggplot(aes(value, group = id, colour = as.factor(mi))) + 
  # geom_density(alpha = .1) + 
  geom_line(aes(alpha = as.factor(mi)), stat = 'density') + 
  geom_line(aes(y = theory, group = NULL), data = theor,
            colour = 'magenta', lwd = 1.1, ) + 
  facet_wrap(~name, scales = 'free', labeller = label_parsed, ncol = 2) +
  labs(title = "Gaussian", y = 'Density', 
       x = '"True" value of random effects',
       colour = bquote(m[i]))  +
  scale_alpha_manual('dummy',
                     values = c(.05, .05, .05), guide = 'none') + 
  scale_colour_manual(values = c(alpha('red', 1), 
                                 alpha('blue', 1),
                                 alpha('darkgreen', 1)))


# Poisson -----------------------------------------------------------------
DP <- diag(c(.25, .05))
bs5 <- lapply(fitsets$`n = 250, mi = 5, fam = poisson`, get.post.b)
bs10 <- lapply(fitsets$`n = 250, mi = 10, fam = poisson`, get.post.b)
bs15 <- lapply(fitsets$`n = 250, mi = 15, fam = poisson`, get.post.b)

allbs5 <- as.data.frame(do.call(rbind, bs5))
allbs10 <- as.data.frame(do.call(rbind, bs10))
allbs15 <- as.data.frame(do.call(rbind, bs15))
allbs <- rbind(allbs5, allbs10, allbs15)
allbs$id <- rep(1:300, each = 250)
allbs$mi <- rep(c(5, 10, 15), each = (250 * 100))

theor <- do.call(rbind,sapply(1:2, function(x){
  b <- allbs[,x,drop=T]
  den <- density(b)
  value <- seq(min(den$x), max(den$x), length.out = 1e3)
  theory <- dnorm(value, mean = 0, sd = sqrt(DP[x,x]))
  out <- data.frame(value = value, theory = theory)
  names(out)[1] <- colnames(allbs)[x]
  pivot_longer(out, -theory)
}, simplify = F))

allbs %>% 
  pivot_longer(`b[0]`:`b[1]`) %>%
  # filter(id %in% c(1:10, 100:110, 200:210)) %>%
  ggplot(aes(value, group = id, colour = as.factor(mi))) + 
  # geom_density(alpha = .1) + 
  geom_line(aes(alpha = as.factor(mi)), stat = 'density') + 
  geom_line(aes(y = theory, group = NULL), data = theor,
            colour = 'magenta', lwd = 1.1, ) + 
  # facet_wrap(~name, scales = 'free', labeller = label_parsed, ncol = 2) +
  facet_grid(name~mi,scales='free')+
  labs(title = "Poisson", y = 'Density', 
       x = '"True" value of random effects',
       colour = bquote(m[i]))  +
  scale_alpha_manual('dummy',
                     values = c(.05, .05, .05), guide = 'none') + 
  scale_colour_manual(values = c(alpha('red', 1), 
                                 alpha('blue', 1),
                                 alpha('darkgreen', 1)))



allbs %>% 
  pivot_longer(`b[0]`:`b[1]`) %>%
  mutate(theor = ifelse(grepl("1", name), dnorm(value, 0, sqrt(.05)),
                        dnorm(value, 0, sqrt(.25)))) %>% 
  arrange(mi) %>%
  mutate(mi_lab = paste0("m[i]==",mi)) %>%
  mutate(mi_lab = factor(mi_lab, levels = c('m[i]==5', 'm[i]==10', 'm[i]==15'))) %>%
  # filter(id %in% c(1:10, 100:110, 200:210)) %>%
  ggplot(aes(value, group = id, colour = mi_lab)) + 
  # geom_density(alpha = .1) + 
  geom_line(aes(alpha = as.factor(mi)), stat = 'density') + 
  geom_line(aes(y = theor, group = NULL),# data = theor,
            colour = 'magenta', lwd = 1.1) + 
  # facet_wrap(~name, scales = 'free', labeller = label_parsed, ncol = 2) +
  facet_grid(name~mi_lab,scales='free', labeller = label_parsed)+
  labs(title = "Poisson", y = 'Density', 
       x = '"True" value of random effects',
       colour = bquote(m[i]))  +
  scale_alpha_manual('dummy',
                     values = c(.05, .05, .05), guide = 'none') + 
  scale_colour_manual(values = c(alpha('red', 1), 
                                 alpha('blue', 1),
                                 alpha('darkgreen', 1)),
                      guide = 'none')

ggsave("~/Downloads/PoissonMVN.png", width = 180, height = 90, units = 'mm')


# Binomial ----------------------------------------------------------------
DB <- matrix(0.25,1,1)
bs5 <- lapply(fitsets$`n = 250, mi = 5, fam = binomial`, get.post.b)
bs10 <- lapply(fitsets$`n = 250, mi = 10, fam = binomial`, get.post.b)
bs15 <- lapply(fitsets$`n = 250, mi = 15, fam = binomial`, get.post.b)

allbs5 <- as.data.frame(do.call(rbind, bs5))
allbs10 <- as.data.frame(do.call(rbind, bs10))
allbs15 <- as.data.frame(do.call(rbind, bs15))
allbs <- rbind(allbs5, allbs10, allbs15)
allbs$id <- rep(1:300, each = 250)
allbs$mi <- rep(c(5, 10, 15), each = (250 * 100))

allbslong <- allbs %>% pivot_longer(`b[0]`)

theor <- do.call(c, with(allbslong, tapply(value, mi, function(x)
  dnorm(x, mean = 0, sd = sqrt(.25)))))
allbslong$theor <- theor

allbslong %>% 
  arrange(mi) %>%
  mutate(mi_lab = paste0("m[i]==",mi)) %>%
  mutate(mi_lab = factor(mi_lab, levels = c('m[i]==5', 'm[i]==10', 'm[i]==15'))) %>%
  # filter(id %in% c(1:10, 100:110, 200:210)) %>%
  ggplot(aes(value, group = id, colour = mi_lab)) + 
  # geom_density(alpha = .1) + 
  geom_line(aes(alpha = as.factor(mi)), stat = 'density') + 
  geom_line(aes(y = theor, group = NULL),# data = theor,
            colour = 'magenta', lwd = 1.1) + 
  # facet_wrap(~name, scales = 'free', labeller = label_parsed, ncol = 2) +
  facet_grid(name~mi_lab,scales='free', labeller = label_parsed)+
  labs(title = "Binomial", y = 'Density', 
       x = '"True" value of random effects',
       colour = bquote(m[i]))  +
  scale_alpha_manual('dummy',
                     values = c(.05, .05, .05), guide = 'none') + 
  scale_colour_manual(values = c(alpha('red', 1), 
                                 alpha('blue', 1),
                                 alpha('darkgreen', 1)),
                      guide = 'none')



