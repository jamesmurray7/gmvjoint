# Restate true parameter values
D.true <- diag(c(0.25, 0.09, 0.30, 0.06,       # Gaussian (1, 2)
            0.20, 0.04, 0.50, 0.09,       # Count    (1, 2)
            2.00))                        # Binary   (1)
ints <- expand.grid(c(1,3,5,7,9), c(1,3,5,7,9))
ints <- ints[ints$Var1!=ints$Var2,]
D.true[cbind(ints$Var1, ints$Var2)] <- 0.125
beta.true <- c(2, -0.1,  0.1, -0.2, 
              -2, 0.1 , -0.1,  0.2,                                         
               2, -0.1,  0.1, -0.2,
               2, -0.1,  0.1, -0.2,
               1,   -1,    1,   -1)
sigma.true <- c(0.16, 0.16)
gamma.true <- c(0.25, -0.25, 0.25, -0.25, 0.30)
zeta.true <- -0.2

target <- c(vech(D.true), beta.true, sigma.true, gamma.true, zeta.true)
target.mat <- t(apply(t(matrix(target)),2,rep,500))
# Parsing fits.
qz <- qnorm(.975)
# Function to extract estimates
extract.estimates <- function(L) sapply(L, function(x){
  setNames(c(vech(x$coeffs$D), x$coeffs$beta, unlist(x$coeffs$sigma)[1:2], x$coeffs$gamma, x$coeffs$zeta),
           names(x$SE))
})
# Function to extract SEs
extract.SE <- function(L) sapply(L, function(x){
  x$SE
})

ests <- extract.estimates(fits)
SEs <- extract.SE(fits)
times <- sapply(fits, function(y) y$elapsed)

# Empirical Means, SDs and average estimated SE.
emp.mean <- rowMeans(ests)
emp.sd <- apply(ests, 1, sd)
avg.se <- rowMeans(SEs)

# 95% coverage probability
lb <- ests - qz * SEs; ub <- ests + qz * SEs
CPs <- rowSums(lb <= target.mat & ub >= target.mat)/500

# MSE / bias
MSEs <- rowMeans((target.mat - ests)^2)
Biases <- rowMeans(ests-target.mat)

# Elapsed time
Elapseds <- times[1,] + times[2,]
Totals <- times[3,]
qn <- quantile(Elapseds)
qn2 <- quantile(Totals)

# Plotting survival parameters --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
source('zzz/theme_csda.R')

 
ests.df <- as.data.frame(t(ests)) %>% 
  select(`gamma_1`:`zeta_bin`) %>% 
  pivot_longer(`gamma_1`:`zeta_bin`) %>% 
  mutate_at("name", ~ ifelse(grepl("gamma", .x), paste0(gsub("_",'[',.x),']'), "zeta")) %>% 
  mutate(target = case_when(
    grepl("1|3", name) ~ 0.25,
    grepl("2|4", name) ~ -0.25,
    grepl("5", name) ~ 0.30,
    T ~ -0.2
  ))

# Option 1, paneled ----
ests.df %>% 
  ggplot(aes(x = name, y = value)) + 
  geom_boxplot(fill = .nice.orange, outlier.alpha = .33, outlier.size = .50, lwd = .25, fatten = 2) + 
  geom_point(aes(y = target), pch = "*", size = 4.95, colour = '#4c98fe') + 
  facet_wrap(~name, scales = 'free', labeller = label_parsed) + 
  labs(y = 'Estimate', x = "")+
  theme_csda() + 
  theme(
    strip.text = element_text(size=9.5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )  

# Option 2, not paneled ----
# ests.df %>% 
#   ggplot(aes(x = name, y = value)) + 
#   geom_boxplot(fill = .nice.orange, outlier.alpha = .33, outlier.size = .50) + 
#   geom_point(aes(y = target), pch = "*", size = 6.5, colour = 'blue') + 
#   labs(y = 'Estimate', x = 'Parameter') +  
#   scale_x_discrete(labels = parse(text=unique(ests.df$name))) + 
#   theme_csda() + 
#   theme(
#     axis.text.x = element_text(size=9.5),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_text(size = 9.5)
#   )

ggsave("~/Downloads/sim2fig.png", width = 140, height = 90, units = 'mm')

tiff(filename = "~/Downloads/sim2fig.tiff", width = 140, height = 90, units = 'mm',
     res=1e3, compression = 'lzw')
last_plot()
dev.off()



# Tabulating for Supplement -----------------------------------------------
D.to.report <- which(vech(D.true) != 0)
to.report <- c(D.to.report, (max(D.to.report)+1):nrow(ests))

ests.to.report <- ests[to.report,]

# Obtain CPs, biases and MSEs.
Emp.Mean <- emp.mean[to.report]; Emp.SE <- emp.sd[to.report]
Mean.SE <- avg.se[to.report]; CP <- CPs[to.report]; Bias <- Biases[to.report]; MSE <- MSEs[to.report]

to.round <- cbind.data.frame(`Empirical Mean` = Emp.Mean, `Empirical SE` = Emp.SE,
                  `Mean SE` = Mean.SE, Bias = Bias, MSE = MSE, CP = CP)
to.round <- apply(to.round, 2, function(x) format(round(x, 3), nsmall = 3))

tab <- cbind.data.frame(parameter = names(CP), to.round, stringsAsFactors=FALSE)
row.names(tab) <- NULL
tab

target.df <- data.frame(
  parameter = unique(tab$parameter),
  true = target[target!=0]
)

tab <- tab %>% 
  left_join(., target.df, "parameter") %>% 
  mutate(parameter2 = case_when(
    grepl("^D\\[", parameter) ~ stringr::str_replace_all(parameter, "\\[", '_{') %>% 
      stringr::str_replace_all(., "\\]", '}') %>% gsub('\\,','',.),
    parameter == "zeta_bin" ~ 'zeta',
    grepl("^Y\\.\\d", parameter) ~ gsub("^Y\\.", 'beta_{', parameter),
    T ~ parameter
  ))
tab$parameter2 <- gsub("\\_\\(Intercept\\)", "0}", tab$parameter2)
tab$parameter2 <- gsub("\\_time", "1}", tab$parameter2)
tab$parameter2 <- gsub("\\_cont", "2}", tab$parameter2)
tab$parameter2 <- gsub("\\_bin", "3}", tab$parameter2)
tab$Parameter <- paste0("$\\", tab$parameter2, "=", format(round(tab$true, 3), nsmall = 3), "$")
tab$`Emp. Mean (SD)` <- paste0(tab$`Empirical Mean`, ' (', tab$`Empirical SE`, ')')

output.tab <- tab %>% 
  select(Parameter, `Emp. Mean (SD)`, `Mean SE`, Bias, MSE, CP)

qn <- round(qn, 3); qn2 <- round(qn2,3)
caption <- paste0("Parameter estimates for five-variate simulation scenario.",
                  " `Emp. Mean (SD)' denotes the average estimated value with the standard deviation of parameter estimates",
                  " and Mean SE the mean standard error calculated for each parameter from each model fit.",
                  " Coverage probabilities are calculated from $\\hbO\\pm1.96\\mathrm{SE}(\\hbO)$. The",
                  " median [IQR] elapsed time for the approximate EM algorithm to converge and standard",
                  " errors calculated was ", qn[3], " [", qn[2], ", ", qn[4], "] seconds.",
                  " Total computation time (e.g. including time taken to obtain initial estimates etc.)",
                  " was ", qn2[3], " [", qn2[2], ", ", qn2[4], "] seconds.")

xt <- xtable::xtable(output.tab, caption = caption,
                     align = c('l', rep("r", ncol(output.tab))))

# Separate-out D[x,x] into appendices manually?
print(xt,
      include.rownames = FALSE,
      sanitize.text.function = identity,
      sanitize.colnames.function = function(x) gsub("\\_\\d?\\d", "", x))

