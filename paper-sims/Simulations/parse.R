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
CPs <- Map(function(estimate, se){
  lb <- estimate - qz * se; ub <- estimate + qz * se
  rowSums(lb <= target.mat & ub >= target.mat)/100
}, estimate = ests, se = SEs)



# Plotting survival parameters --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
source('zzz/theme_csda.R')

all.ests <- do.call(rbind, lapply(seq_along(ests), function(x){
  current <- names(ests)[x]
  mi <- stringr::str_extract(current, 'mi\\s\\=\\s\\d?\\d') %>% 
        gsub("mi\\s\\=\\s",'',.) %>% as.integer
  rate <- stringr::str_extract(current, "\\d\\d\\%$") %>% 
        gsub("\\%", "", .) %>% as.integer
  
  survs <- as.data.frame(t(ests[[x]][(nrow(ests[[x]])-3):nrow(ests[[x]]),]))
  info <- data.frame(r = rep(mi, nrow(survs)),
                     omega = rep(rate/100, nrow(survs)))
  cbind(info, survs)
}))

all.ests %>% 
  pivot_longer(`gamma_1`:`zeta_bin`) %>% 
  mutate(name = case_when(
    name == "gamma_1" ~ "gamma[1]",
    name == "gamma_2" ~ "gamma[2]",
    name == "gamma_3" ~ "gamma[3]",
    T ~ "zeta"
  )) %>% 
  mutate(target = case_when(
    name == "gamma[1]" ~ .5,
    name == "gamma[2]" ~ -.5,
    name == "gamma[3]" ~ .5,
    T ~ -.2
  )) %>% 
  mutate(r = factor(r, levels = c(5,10,15)),
         omega = factor(omega, levels = c(.1, .3, .5))) %>% 
  ggplot(aes(x = omega, y = value, fill = r)) + 
  geom_hline(aes(yintercept = target), lty = 5, colour = 'grey5', alpha = .5) + 
  geom_boxplot(outlier.alpha = .33) + 
  facet_wrap(~name, scales = 'free_y', labeller = label_parsed) + 
  labs(fill = "Maximal profile length", 
       x = "Failure rate",
       y = "Estimate") + 
  theme_csda() + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme(strip.text = element_text(size=9.5),
        legend.position = 'bottom')

ggsave("~/Downloads/sim1fig.png", width = 140, height = 90, units = 'mm')

tiff(filename = "~/Downloads/sim1fig.tiff", width = 140, height = 90, units = 'mm',
     res=1e3, compression = 'lzw')
last_plot()
dev.off()
  