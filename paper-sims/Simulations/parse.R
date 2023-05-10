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
times <- lapply(fits, function(x) sapply(x, function(y) y$elapsed))

# Empirical Means, SDs and average estimated SE.
emp.mean <- lapply(ests, rowMeans)
emp.sd <- lapply(ests, function(x) apply(x, 1, sd))
avg.se <- lapply(SEs, rowMeans)

# 95% coverage probability
CPs <- Map(function(estimate, se){
  lb <- estimate - qz * se; ub <- estimate + qz * se
  rowSums(lb <= target.mat & ub >= target.mat)/100
}, estimate = ests, se = SEs)

MSEs <- Map(function(estimate){
  rowMeans((target.mat - estimate)^2)
}, estimate = ests)

Biases <- Map(function(estimate){
  rowMeans(estimate-target.mat)
}, estimate = ests)

Elapseds <- Map(function(TIME){
  EM <- TIME[1,]
  PP <- TIME[2,]
  EM + PP
}, TIME = times)

Totals <- Map(function(TIME){
  TIME[3,]
}, TIME = times)

iter.per.second <- Map(function(TIME){
  TIME[4,]/TIME[1,]
}, TIME = times)

df.Elapsed <- expand.grid(r = c(5, 10, 15),
                          omega = c(0.1, 0.3, 0.5))
df.Elapsed <- cbind(df.Elapsed, t(apply(df.Elapsed,1,function(x){
  r <- x[1]; fail <- paste0(100*x[2], "%")
  .lookup <- paste0("n = 250, mi = ", r, ", failure = ", fail)
  qn <- quantile(Elapseds[[.lookup]])
  c(qn[2], qn[3], qn[4])
})))

df.Totals <- expand.grid(r = c(5, 10, 15),
                          omega = c(0.1, 0.3, 0.5))
df.Totals <- cbind(df.Totals, t(apply(df.Totals,1,function(x){
  r <- x[1]; fail <- paste0(100*x[2], "%")
  .lookup <- paste0("n = 250, mi = ", r, ", failure = ", fail)
  qn <- quantile(Totals[[.lookup]])
  c(qn[2], qn[3], qn[4])
})))

sapply(iter.per.second, mean)

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
  geom_boxplot(outlier.alpha = .33,
               outlier.size = .50) + 
  facet_wrap(~name, scales = 'free_y', labeller = label_parsed) + 
  labs(fill = expression(r~"="),# "r =",#Maximal profile length", 
       x = expression(omega),
       y = "Estimate") + 
  theme_csda() + 
  scale_fill_brewer(palette = 'YlOrRd') + 
  theme(strip.text = element_text(size=9.5),
        legend.title = element_text(size=8),
        legend.position = 'bottom')

ggsave("~/Downloads/sim1fig.png", width = 140, height = 90, units = 'mm')

tiff(filename = "~/Downloads/sim1fig.tiff", width = 140, height = 90, units = 'mm',
     res=1e3, compression = 'lzw')
last_plot()
dev.off()



# Tabulating for Supplement -----------------------------------------------
all.ests <- do.call(rbind, lapply(seq_along(ests), function(x){
  current <- names(ests)[x]
  mi <- stringr::str_extract(current, 'mi\\s\\=\\s\\d?\\d') %>% 
    gsub("mi\\s\\=\\s",'',.) %>% as.integer
  rate <- stringr::str_extract(current, "\\d\\d\\%$") %>% 
    gsub("\\%", "", .) %>% as.integer
  
  tests <- t(ests[[x]])
  tests <- tests[c(1,3,5,6,10,12,13,15:ncol(tests)), c(1,3,5,6,10,12,13,15:ncol(tests))]
  info <- data.frame(r = rep(mi, nrow(tests)),
                     omega = rep(rate/100, nrow(tests)))
  cbind(info, tests)
}))
  
means <- all.ests %>% 
  group_by(omega, r) %>% 
  mutate(across(`D[1,1]`:`zeta_bin`, ~mean(.x))) %>% 
  distinct() %>% 
  ungroup

sds <- all.ests %>% 
  group_by(omega, r) %>% 
  mutate(across(`D[1,1]`:`zeta_bin`, ~sd(.x))) %>% 
  distinct() %>% 
  ungroup

msd <- left_join(
  means %>% 
    pivot_longer(cols = `D[1,1]`:`zeta_bin`,
                 names_to = "parameter", values_to = "mean"),
  sds %>% 
    pivot_longer(cols = `D[1,1]`:`zeta_bin`,
                 names_to = "parameter", values_to = "sd"),
  by = c('r', 'omega', 'parameter')
) %>% 
  mutate_at(vars(mean, sd), ~ format(round(.x, 3), nsmall = 3)) %>% 
  mutate(`Empirical mean (SE)` = paste0(mean, ' (', sd, ')'))

# Obtain CPs, biases and MSEs.
all.tab.ests <- do.call(rbind, lapply(seq_along(CPs), function(x){
  current <- names(CPs)[x]
  mi <- stringr::str_extract(current, 'mi\\s\\=\\s\\d?\\d') %>% 
    gsub("mi\\s\\=\\s",'',.) %>% as.integer
  rate <- stringr::str_extract(current, "\\d\\d\\%$") %>% 
    gsub("\\%", "", .) %>% as.integer
  
  Emp.Mean <- emp.mean[[x]]; Emp.SE <- emp.sd[[x]]
  Mean.SE <- avg.se[[x]]; CP <- CPs[[x]]; Bias <- Biases[[x]]; MSE <- MSEs[[x]]
  Emp.Mean <- Emp.Mean[c(1,3,5,6,10,12,13,15:length(Emp.Mean))]
  Emp.SE <- Emp.SE[c(1,3,5,6,10,12,13,15:length(Emp.SE))]
  Mean.SE <- Mean.SE[c(1,3,5,6,10,12,13,15:length(Mean.SE))]
  CP <- CP[c(1,3,5,6,10,12,13,15:length(CP))]
  Bias <- Bias[c(1,3,5,6,10,12,13,15:length(Bias))]
  MSE <- MSE[c(1,3,5,6,10,12,13,15:length(MSE))]
  
  to.round <- cbind.data.frame(`Empirical Mean` = Emp.Mean, `Empirical SE` = Emp.SE,
                    `Mean SE` = Mean.SE, Bias = Bias, MSE = MSE)
  to.round <- apply(to.round, 2, function(x) format(round(x, 3), nsmall = 3))
  
  info <- data.frame(r = rep(mi, length(CP)),
                     omega = rep(rate/100, length(CP)))
  out <- cbind.data.frame(info, parameter = names(CP), to.round, CP = CP, stringsAsFactors=FALSE)
  row.names(out) <- NULL
  out
}))

low <- all.tab.ests[all.tab.ests$omega==.1,]
med <- all.tab.ests[all.tab.ests$omega==.3,]
high <- all.tab.ests[all.tab.ests$omega==.5,]

target.df <- data.frame(
  parameter = unique(low$parameter),
  true = target[target!=0]
)

df.to.xtab <- function(tab){
  # Clean parameter names
  
  med.iqr.elapsed <- df.Elapsed[df.Elapsed$omega == tab$omega[1], ]
  med.iqr.total <- df.Totals[df.Totals$omega == tab$omega[1],]
  med.iqr.report <- apply(med.iqr.elapsed,1,function(x){
    r <- x[1]
    med <- format(round(x[4], 3), nsmall = 3)
    p25 <- format(round(x[3], 3), nsmall = 3)
    p75 <- format(round(x[5], 3), nsmall = 3)
    
    paste0(med, ' [', p25, ', ', p75, "] seconds for $r=", r, "$")
  })
  med.iqr.total.report <- apply(med.iqr.total,1,function(x){
    r <- x[1]
    med <- format(round(x[4], 3), nsmall = 3)
    p25 <- format(round(x[3], 3), nsmall = 3)
    p75 <- format(round(x[5], 3), nsmall = 3)
    
    paste0(med, ' [', p25, ', ', p75, "] seconds for $r=", r, "$")
  })
  
  tab <- tab %>% 
    select(-omega) %>% 
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
    select(r, Parameter, `Emp. Mean (SD)`, `Mean SE`, Bias, MSE, CP)
  
  output.tab <- output.tab %>% 
    pivot_wider(id_cols = Parameter,
                names_from = r, names_vary='slowest',
                values_from = `Emp. Mean (SD)`:CP)
  
  caption <- paste0("Parameter estimates for $\\omega=", med.iqr.elapsed$omega[1], "$",
                    " for differing maximal longitudinal profile lengths $r$. `Emp. Mean (SD)'",
                    " denotes the average estimated value with the standard deviation of parameter estimates",
                    " and Mean SE the mean standard error calculated for each parameter from each model fit.",
                    " Coverage probabilities are calculated from $\\hbO\\pm1.96\\mathrm{SE}(\\hbO)$. The",
                    " median [IQR] elapsed times for the approximate EM algorithm to converge and standard",
                    " errors calculated was ", med.iqr.report[1], ", ", med.iqr.report[2], " and ", med.iqr.report[3], ".",
                    " \\textit{Total} computation time was ", med.iqr.total.report[1], ", ", med.iqr.total.report[2], " and ", med.iqr.total.report[3], ".")
  
  xt <- xtable::xtable(output.tab, caption = caption,
                       align = c('l', rep("r", ncol(output.tab))))
  
  addtorow <- list()
  addtorow$pos <- list(-1)
  addtorow$command <- c("& \\multicolumn{5}{c}{$r=5$} & \\multicolumn{5}{c}{$r=10$} & \\multicolumn{5}{c}{$r=15$}\\\\")
  
  print(xt,
        include.rownames = FALSE,
        add.to.row = addtorow,
        sanitize.text.function = identity,
        sanitize.colnames.function = function(x) gsub("\\_\\d?\\d", "", x))
}

df.to.xtab(low)
df.to.xtab(med)
df.to.xtab(high)
