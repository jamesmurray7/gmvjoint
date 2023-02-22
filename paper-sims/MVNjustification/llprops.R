.sim <- function(ntms, theta, family){
  a <- simData(n = 250, ntms = ntms, theta = theta,
          beta = t(c(2, -0.1, 0.1, -0.2)),
          D =  diag(c(0.25, 0.05)),
          zeta = c(0, -0.2),
          family = as.list(family),
          gamma = 0.5,
          return.ranefs = TRUE)
  list(a$data, a$ranefs)
}

# Some good thetas to use.
theta70ish <- c(-1.5, 0.1)
theta50ish <- c(-2, 0.1)
theta20ish <- c(-3, 0.1)

# Candidate m_i, theta class and family
to.sim <- expand.grid(mi = c(5, 10, 15),
                      theta = c('low', 'medium', 'high'),
                      family = c('gaussian', 'poisson', 'binomial'))

# Simulate
sim.sets <- apply(to.sim, 1, function(x){
  m <- as.numeric(x[1]); cl <- x[2]; fam <- x[3]
  if(cl == 'low'){
    theta <- theta20ish
  }else if(cl == 'medium'){
    theta <- theta50ish
  }else{
    theta <- theta70ish
  }
  
  replicate(100, .sim(m, theta, fam), simplify = F)
})
names(sim.sets) <- sapply(apply(to.sim, 1, paste0, collapse=','), trimws)

# Function to get log-likelihoods
getL <- function(x, family){
  data <- x[[1]]; b <- x[[2]]
  sdata <- data[!duplicated(data$id),]
  X <- model.matrix(~time+cont+bin, data)
  Z <- model.matrix(~time, data)
  Y <- data$Y.1
  eta <- X %*% c(2, -0.1, 0.1, -0.2) + rowSums(Z * b[data$id, ])
  longL <- switch(family, 
              gaussian = dnorm(Y, eta, sd = sqrt(.16), log = T),
              poisson = dpois(Y, exp(eta), log = T),
              binomial = dbinom(Y, 1, plogis(eta), log = T))
  bL <- mvtnorm::dmvnorm(b, sigma = diag(c(0.25, 0.05)), log = T)
  survL <- c(logLik(coxph(Surv(survtime, status) ~ bin, sdata)))
  longL <- sum(longL)
  bL <- sum(bL)
  list(long = longL, ranef = bL, surv = survL, 
       long.res = longL + bL,
       l = survL + longL + bL)
}

Lsets <- setNames(vector('list', length(sim.sets)), names(sim.sets))
fams <- as.character(to.sim$family)
for(p in 1:length(sim.sets)){
  Lsets[[p]] <- lapply(sim.sets[[p]], getL, fams[p])
  cat(sprintf("%d/%d\r", p, length(sim.sets)))
}

# Some neat way of plotting?
props <- lapply(Lsets, function(y)
  do.call(rbind, lapply(y, function(x) c(long = x$long.res, surv = x$surv, total = x$l))))

props.df <- as.data.frame(do.call(rbind, props))
props.df$name <- rep(names(sim.sets), each=100)
props.df$longprop <- props.df$long/props.df$total
props.df$survprop <- props.df$surv/props.df$total

boxplot(props.df$longprop ~ props.df$name)

library(ggplot2)
library(dplyr)
library(tidyr)

props.df %>% 
  mutate(failure = gsub(',','',stringr::str_extract(pattern = '\\,\\w+\\,', name)),
         family = gsub(',','',stringr::str_extract(pattern = '\\w+$', name)),
         mi = stringr::str_extract(name, '\\d+')) %>%
  mutate(mi = factor(mi, levels = c('5','10','15')),
         failure = factor(failure, levels = c('low', 'medium', 'high'))) %>%
  select(failure, family, mi, longprop) %>% 
  pivot_longer(longprop) %>%
  ggplot(aes(x = mi, y = value)) + 
  geom_hline(aes(yintercept=0.5), colour = 'magenta', lty = 3) + 
  geom_boxplot() + 
  facet_grid(failure~family, scales = 'free_y') + 
  theme_light() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour='grey3'),
    panel.grid.major.x = element_blank()
  ) + 
  xlab(expression(m[i])) + 
  ylab(expression("\u2113"["long"]/"\u2113"["total"]))

         