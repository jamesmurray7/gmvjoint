rm(list=ls())
bigD <- TRUE
# Data --------------------------------------------------------------------
data.dir <- '/data/c0061461/fits_ghnodes/'
if(bigD) data.dir <- gsub("\\/$", "_largeD/", data.dir)
RDs <- dir(data.dir, pattern = '.RD')
RDs <- RDs[-which(RDs%in%c("datas.RData", "tab.txt"))]

# Targets -----------------------------------------------------------------
if(!bigD) Dtr <- diag(c(0.25, 0.09, 0.50, 0.10, 2.00)) else Dtr <- diag(c(2, .4, 1.5, 0.33, 3.00))

inds <- expand.grid(c(1,3,5), c(1,3,5))
inds <- inds[inds$Var1!=inds$Var2,]
Dtr[cbind(inds$Var1, inds$Var2)] <- if(!bigD) 0.25 else 0.50
nmD <- paste0('D[',apply(which(lower.tri(Dtr, T), arr.ind= T),1,paste,collapse=','),']')

# Parameters 
betatr <- rbind(                         # Fixed effects
  c(2, -0.1, 0.1, -0.2),
  c(2, -0.1, 0.1, -0.2),
  c(1, -1, 1, -1)
)
nmb <- paste0(rep(paste0("Y.", 1:3), each = 4), c("_(Intercept)", "_time", "_cont", "_bin"))

targets <- setNames(c(vech(Dtr),            #vech(D)
                      c(apply(betatr,1,c)), #beta
                      0.16,                 #sigma2 
                      0.50, -0.5, 0.5,      #gammas
                     -0.20),
                    c(nmD, nmb,
                      "sigma^2_1", paste0("gamma_",1:3), "zeta_bin"))

target.mat <- matrix(targets, nrow = length(targets), ncol = 100)
row.names(target.mat) <- names(targets)

# Get _all_ estimates -----------------------------------------------------
# (This is for __tabulation__)
qz <- qnorm(.975)
getEsts <- function(X){
  # X a list of lists
  num.null <- sum(sapply(X, is.null))
  SEs <- sapply(X, function(x){
    if(is.null(x)) 
      return(NULL)
    else
      return(x$SE)
  })
  if(num.null > 0) SEs <- do.call(cbind, SEs)
  nm <- row.names(SEs)
  ests <- sapply(X, function(x){
    if(is.null(x))
      return(NULL)
    else
      return(setNames(c(vech(x$coeffs$D), x$coeffs$beta, x$coeffs$sigma[[1]], x$coeff$gamma, x$coeff$zeta),
                      nm))
  })
  if(num.null > 0) ests <- do.call(cbind, ests)
  lbs <- ests - qz * SEs
  ubs <- ests + qz * SEs
  
  # Empirical mean (SD), average SE.
  Emp.Mean <- rowMeans(ests)
  Emp.SD <- apply(ests, 1, sd)
  Avg.SE <- rowMeans(SEs)
  # Bias
  TM <- target.mat[,1:ncol(ests)]
  Bias <- rowMeans(ests - TM)
  MSE <- rowMeans((TM - ests)^2)
  # CP
  CP <- rowSums(lbs <= TM & ubs >= TM)/ncol(ests)
  list(
    Emp.Mean = Emp.Mean,
    Emp.SD = Emp.SD,
    Avg.SE = Avg.SE,
    Bias = Bias,
    MSE = MSE,
    CP = CP
  )
}


# Sequence along files ----------------------------------------------------
# and tabulate
.f <- function(x, d = 3) format(round(x, d), nsmall = d)
tabs <- lapply(seq_along(RDs), function(i){
  f <- RDs[i]
  file <- paste0(data.dir,f)
  cat(sprintf("Current file: %s\n", file))
  nodes <- as.numeric(gsub('\\.RData', '',gsub('n\\_', '', f)))
  assign("x", get(load(file)))
  fits <- x#[[1]]
  num.null <- sum(sapply(fits, is.null))
  cat(sprintf("%d NULL fits.\n", num.null))
  # Get "S"ummary
  S <- getEsts(fits)
  MSD <- paste0(.f(S$Emp.Mean), ' (', .f(S$Emp.SD), ')')
  ptar <- paste0(names(targets), ' = ', .f(targets))
  out <- data.frame(nodes = as.integer(.f(nodes, 0)), parameter = ptar, `MeanSD` = MSD, 
                    SE = .f(S$Avg.SE), Bias = .f(S$Bias), MSE = .f(S$MSE),
                    CP = .f(S$CP, 2), stringsAsFactors = F)
  row.names(out) <- NULL
  cat(cli::bg_cyan(cli::col_red('Done!')))
  cat("\n\n")
  out
})

tabs2 <- do.call(rbind,tabs)
tabs2 <- dplyr::arrange(tabs2, nodes) %>% 
  filter(grepl("Y\\.3|sigma|gamma|zeta", parameter)) # Only these ones are taken with quadrature.

tabs2


# WIDE version
tabsw <- tidyr::pivot_wider(tabs2, 
                   id_cols = parameter,
                   names_from = nodes, values_from = MeanSD:CP,
                   names_vary = 'slowest')
# Clean parameters up.
library(dplyr)
tabsw <- tabsw %>% 
  mutate(parameter2 = case_when(
    grepl("^D\\[", parameter) ~ stringr::str_replace_all(parameter, "\\[", '_{') %>% 
      stringr::str_replace_all(., "\\]", '}') %>% gsub('\\,','',.),
    grepl("zeta_bin", parameter) ~ 'zeta = -0.200',
    grepl("^Y\\.\\d", parameter) ~ gsub("^Y\\.", 'beta_{', parameter),
    T ~ parameter
  ))
tabsw$parameter2 <- gsub("\\_\\(Intercept\\)", "0}", tabsw$parameter2)
tabsw$parameter2 <- gsub("\\_time", "1}", tabsw$parameter2)
tabsw$parameter2 <- gsub("\\_cont", "2}", tabsw$parameter2)
tabsw$parameter2 <- gsub("\\_bin", "3}", tabsw$parameter2)
tabsw$Parameter <- paste0("$\\", tabsw$parameter2, '$')
tabsw$Parameter <- gsub("\\_1","",tabsw$Parameter)
tabsw$Parameter <- gsub("=\\s\\s", "= ", tabsw$Parameter)

nms <- names(tabsw)
nodes <- unique(tabs2$node)

library(xtable)
make.xt <- function(x, comm = NULL, ...){
  xt <- xtable(x)
  if(!is.null(comm)){
    addtorow <- list()
    addtorow$pos <- list(-1)
    addtorow$command <- comm  
  }
  print(xt,
        include.rownames = FALSE,
        add.to.row = if(!is.null(comm)) addtorow else NULL,
        sanitize.text.function = identity,...)
  # sanitize.colnames.function = function(y) gsub('\\.', ' ', y))
}

# 3 sets of 5?
# dfs.to.tabulate <- c(2,3,5,10,30,50)
nodes.to.tabulate <- unique(tabs2$nodes)
stats <- c("MeanSD", "MSE", "CP")

cols <- paste0(stats,'_',rep(nodes.to.tabulate, each = length(stats)))

tabsw[,c('Parameter', cols)]

# cols <- list()
# for(j in 1:tab.cols){
#   rowdf <- rows[[j]]
#   out <- vector('list', length(rowdf))
#   for(d in seq_along(rowdf)){
#     cols.to.extract <- which(stringr::str_extract(nms, '\\_\\d?\\d?\\d$') %in% paste0('_', rowdf[d]))
#     this.df <- as.data.frame(cbind(Pcol, tabsw[, cols.to.extract]))
#     names(this.df) <- gsub('\\_.*', '', names(this.df))
#     names(this.df)[which(names(this.df) == 'MeanSD')] <- "Mean (SE)"
#     names(this.df)[which(names(this.df) == 'Pcol')] <- "Parameter"
#     out[[d]] <- this.df
#   }
#   cols[[j]] <- do.call(rbind, out)
# }

make.xt(tabsw[,c('Parameter', cols)])

write.table(make.xt(tab),
            file = paste0(data.dir, 'tab.txt'),
            quote = F, )


# Plotting coverages ------------------------------------------------------
library(ggplot2)

tabsplot <- tabs2 %>%
  select(nodes, parameter, Bias, MSE, CP) %>%
  mutate(parameter = trimws(gsub('\\=.*$', '', parameter))) %>%
  mutate(parameter2 = case_when(
    grepl("^D\\[", parameter) ~ stringr::str_replace_all(parameter, "\\[", '[') %>%
      stringr::str_replace_all(., "\\]", ']') %>% gsub('\\,','',.),
    grepl("zeta_bin", parameter) ~ 'zeta',
    grepl("gamma", parameter) ~ paste0('gamma[', readr::parse_number(parameter), ']'),
    grepl("^Y\\.\\d", parameter) ~ gsub("^Y\\.", 'beta[', parameter),
    parameter == "sigma^2_1" ~ "sigma[epsilon[1]]^2",
    T ~ parameter
  ))
tabsplot$parameter2 <- gsub("\\_\\(Intercept\\)", "0]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_time", "1]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_cont", "2]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_bin", "3]", tabsplot$parameter2)

# Create 

tabsplot %>% 
  mutate_at(c("MSE", "Bias", "CP"), as.numeric) %>%
  ggplot(aes(x = nodes, y = abs(Bias))) + 
  # geom_hline(aes(yintercept=0), lty=3)+
  geom_line()+
  scale_x_continuous(name = expression(rho), 
                     breaks = unique(tabsplot$nodes)) +
  facet_wrap(~parameter2,ncol=4,nrow=2, scales = 'free_y', labeller = label_parsed)+
  theme_light() +
  labs(y = 'Absolute bias') + 
  theme(
    strip.background = element_blank(),
    # strip.text = element_text(colour = "black", size = 10, face = 'bold'),
    strip.text = element_text(colour = "black"),
    legend.position = 'none',
    # axis.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title = element_text(size = 5 * 1.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(paste0(data.dir,'plot.png'),
       device = 'png',
       width = 140, height = 90, units = 'mm')


# These quite boring, what about 2.5/97.5%? -------------------------------

newfn <- function(X){
  #x a list 
  # X a list of lists
  num.null <- sum(sapply(X, is.null))
  SEs <- sapply(X, function(x){
    if(is.null(x)) 
      return(NULL)
    else
      return(x$SE)
  })
  if(num.null > 0) SEs <- do.call(cbind, SEs)
  nm <- row.names(SEs)
  ests <- sapply(X, function(x){
    if(is.null(x))
      return(NULL)
    else
      return(setNames(c(vech(x$coeffs$D), x$coeffs$beta, x$coeffs$sigma[[1]], x$coeff$gamma, x$coeff$zeta),
                      nm))
  })
  if(num.null > 0) ests <- do.call(cbind, ests)
  lbs <- ests - qz * SEs
  ubs <- ests + qz * SEs
  
  # Bias
  TM <- target.mat[,1:ncol(ests)]
  # NB this is abs error, not bias!
  Bias <- abs(TM - ests)                       # 2.5%, 25%, 50%, 75%, 97.5%
  Biassumm <- apply(Bias, 1, quantile, probs = c(.025, .25, .5, .75, .975))
  MSE <- ((TM - ests)^2)
  MSEsumm <- apply(MSE, 1, quantile, probs = c(.025, .25, .5, .75, .975))
  # CP
  CP <- rowSums(lbs <= TM & ubs >= TM)/ncol(ests)
  list(
    Bias = Bias,
    Biassumm = Biassumm,
    MSE = MSE,
    MSEsumm = MSEsumm,
    CP = CP
  )
}

to.plot <- lapply(seq_along(RDs), function(i){
  f <- RDs[i]
  file <- paste0(data.dir,f)
  cat(sprintf("Current file: %s\n", file))
  nodes <- as.numeric(gsub('\\.RData', '',gsub('n\\_', '', f)))
  assign("x", get(load(file)))
  fits <- x#[[1]]
  num.null <- sum(sapply(fits, is.null))
  cat(sprintf("%d NULL fits.\n", num.null))
  # Get "S"ummaries
  S <- newfn(fits)
  
  # bias ----
  Bias <- t(S$Biassumm) %>% as.data.frame() %>% tibble::rownames_to_column("parameter")
  names(Bias)[-1] <- paste0("Bias", names(Bias)[-1])
  # MSE ----
  MSE <- t(S$MSEsumm) %>% as.data.frame() %>% tibble::rownames_to_column("parameter")
  names(MSE)[-1] <- paste0("MSE", names(MSE)[-1])
  # CP ----
  CP <- as.data.frame(t(t(S$CP))) %>% tibble::rownames_to_column("parameter")
  names(CP)[-1] <- "CP"
  
  out <- left_join(Bias, MSE, 'parameter') %>% left_join(., CP, 'parameter')
  out$nodes <- nodes
  
  parameters.to.keep <<- grepl("Y\\.3|sigma|gamma|zeta", out$parameter)
  
  cat(cli::bg_cyan(cli::col_red('Done!')))
  cat("\n\n")
  # out[parameters.to.keep,]
  out
})

to.plot2 <- do.call(rbind, to.plot) %>% arrange(nodes)
# Quad parameters, bias
to.plot2 %>% 
  filter(grepl("Y\\.3|sigma|gamma|zeta", parameter), nodes < 90) %>% 
  mutate(parameter = case_when(
    grepl("gamma", parameter) ~ paste0("gamma[", stringr::str_extract(parameter, '\\d'), ']'),
    grepl("sigma", parameter) ~ "sigma[epsilon]^2",
    parameter == 'zeta_bin'  ~ "zeta",
    T ~ gsub("Y\\.3", "beta[3", gsub("\\_", "", parameter))
  )) %>% 
  mutate(
    parameter = stringr::str_replace(parameter, "\\(Intercept\\)", "0]"),
    parameter = stringr::str_replace(parameter, "time", "1]"),
    parameter = stringr::str_replace(parameter, "cont", "2]"),
    parameter = stringr::str_replace(parameter, "bin", "3]")
  ) %>% 
  ggplot(aes(x = nodes, y = `Bias50%`)) +
  geom_line() + 
  geom_point(pch = 20) + 
  # The quantile lines just make this harder to see, so remove for now...
  # geom_line(aes(y = `Bias25%`), lty = 3, col = 'grey') + 
  # geom_line(aes(y = `Bias75%`), lty = 3, col = 'grey') + 
  facet_wrap(~parameter, scales = 'free_y', labeller = label_parsed) + 
  theme_light() + 
  labs(y = "Median absolute error", x = expression(rho)) + 
  scale_x_continuous(breaks = unique(to.plot2$nodes)) + 
  theme(
    strip.background = element_blank(),
    # strip.text = element_text(colour = "black", size = 10, face = 'bold'),
    strip.text = element_text(colour = "black"),
    legend.position = 'none',
    # axis.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title = element_text(size = 5 * 1.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(paste0(data.dir, "errorplot.png"),
       device = 'png',
       width = 140, height = 90, units = 'mm')
