bigD <- TRUE # change if small
# Data --------------------------------------------------------------------
data.dir <- '/data/c0061461/fits_tdistn/'
if(bigD) data.dir <- gsub('\\/$', '_largeD/', data.dir)
RDs <- dir(data.dir)

# Targets -----------------------------------------------------------------
if(bigD) vD <- c(2,0,.4) else vD <- c(.25, 0, 0.05)
targets <- setNames(c(vD,
                      2.00, 0.33,-0.50, 0.25,
                      0.16,
                      0.50,
                     -0.20),
                    c("D[1,1]", "D[2,1]", "D[2,2]",
                      "Y.1_(Intercept)", "Y.1_time", "Y.1_cont", "Y.1_bin",
                      "sigma^2_1", "gamma_1", "zeta_bin"))

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
  df <- as.numeric(gsub('\\.RData', '',gsub('df\\_', '', f)))
  assign("x", get(load(file)))
  fits <- x[[1]]
  num.null <- sum(sapply(fits, is.null))
  cat(sprintf("%d NULL fits.\n", num.null))
  # Get "S"ummary
  S <- getEsts(fits)
  MSD <- paste0(.f(S$Emp.Mean), ' (', .f(S$Emp.SD), ')')
  ptar <- paste0(names(targets), ' = ', .f(targets))
  out <- data.frame(df = as.integer(.f(df, 0)), parameter = ptar, `MeanSD` = MSD, 
                    SE = .f(S$Avg.SE), Bias = .f(S$Bias), MSE = .f(S$MSE),
                    CP = .f(S$CP, 2), stringsAsFactors = F)
  row.names(out) <- NULL
  cat(cli::bg_cyan(cli::col_red('Done!')))
  cat("\n\n")
  out
})

tabs2 <- do.call(rbind,tabs)
tabs2 <- dplyr::arrange(tabs2, df)


# WIDE version
tabsw <- tidyr::pivot_wider(tabs2, 
                   id_cols = parameter,
                   names_from = df, values_from = MeanSD:CP,
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
dfs <- unique(tabs2$df)

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
dfs.to.tabulate <- c(2,3,5,10,30,50)
tab.rows <- 3
tab.cols <- 2
rows <- list(dfs.to.tabulate[1:3], dfs.to.tabulate[4:length(dfs.to.tabulate)])
Pcol <- tabsw$Parameter                                # so just check and move on!

cols <- list()
for(j in 1:tab.cols){
  rowdf <- rows[[j]]
  out <- vector('list', length(rowdf))
  for(d in seq_along(rowdf)){
    cols.to.extract <- which(stringr::str_extract(nms, '\\_\\d?\\d?\\d$') %in% paste0('_', rowdf[d]))
    this.df <- as.data.frame(cbind(Pcol, tabsw[, cols.to.extract]))
    names(this.df) <- gsub('\\_.*', '', names(this.df))
    names(this.df)[which(names(this.df) == 'MeanSD')] <- "Mean (SE)"
    names(this.df)[which(names(this.df) == 'Pcol')] <- "Parameter"
    out[[d]] <- this.df
  }
  cols[[j]] <- do.call(rbind, out)
}

lhs <- cols[[1]]; rhs <- cols[[2]]
nr.rhs <- nrow(rhs); nr.lhs <- nrow(lhs)
while(nr.rhs < nr.lhs){
  rhs <- suppressWarnings(rbind(rhs, rep("{}", 6)))
  nr.rhs <- nrow(rhs)
}

tab <- cbind(lhs,rhs[,-1])
make.xt(tab)
write.table(make.xt(tab), 
            file = paste0(data.dir, 'tab.txt'),
            quote=F)

# Plotting coverages ------------------------------------------------------
library(ggplot2)

tabsplot <- tabs2 %>% 
  select(df, parameter, CP) %>% 
  mutate(parameter = trimws(gsub('\\=.*$', '', parameter))) %>% 
  mutate(parameter2 = case_when(
    grepl("^D\\[", parameter) ~ stringr::str_replace_all(parameter, "\\[", '[') %>% 
      stringr::str_replace_all(., "\\]", ']') %>% gsub('\\,','',.),
    grepl("zeta_bin", parameter) ~ 'zeta',
    grepl("^Y\\.\\d", parameter) ~ gsub("^Y\\.", 'beta[', parameter),
    parameter == "sigma^2_1" ~ "sigma^2",
    parameter == "gamma_1"  ~ "gamma",
    T ~ parameter
  ))
tabsplot$parameter2 <- gsub("\\_\\(Intercept\\)", "0]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_time", "1]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_cont", "2]", tabsplot$parameter2)
tabsplot$parameter2 <- gsub("\\_bin", "3]", tabsplot$parameter2)

  
tabsplot %>% 
  mutate_at('CP', as.numeric) %>% 
  ggplot(aes(x = df, y = CP)) + 
  geom_hline(aes(yintercept = 0.95), lty = 5, col = 'lightgray') + 
  geom_line(lwd = 1.2) + 
  labs(y = '95% Coverage Probability') + 
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 1.00)) + 
  scale_x_log10(expression(nu)) + 
  facet_wrap(~parameter2, labeller = label_parsed, ncol = 5, nrow = 2)+
  theme_light() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 12, face = 'bold'),
    legend.position = 'none',
    axis.text.y = element_text(size = 6, angle = 45),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(paste0(data.dir, 'CP.png'), device = 'png',
       width = 190, height = 150, units = 'mm')


# RE Comparison -----------------------------------------------------------
# Point estimate vs actual value

ET <- lapply(seq_along(RDs), function(i){
  f <- RDs[i]
  file <- paste0(data.dir,f)
  cat(sprintf("Current file: %s\n", file))
  df <- as.numeric(gsub('\\.RData', '',gsub('df\\_', '', f)))
  assign("x", get(load(file)))
  fits <- x[[1]]
  true.REs <- x[[2]]
  num.null <- sum(sapply(fits, is.null))
  cat(sprintf("%d NULL fits.\n", num.null))
  
  Est.trues <- lapply(seq_len(100), function(x){
    fx <- fits[[x]]
    bx <- true.REs[[x]]
    if(is.null(fx)){ 
      return(NULL)
    }else{
     RE <- fx$REs
     attr(RE, 'Var') <- NULL
     attr(RE, 'vcov') <- NULL
     ids <- 1:nrow(RE)
     est <- data.frame(id = ids, df = df, what = 'Estimate', RE)
     colnames(bx) <- colnames(RE)
     true <- data.frame(id = ids, df = df, what = 'True', bx)
     out <- rbind(est, true)
     out$simnum <- x
     return(out)
    }
  })
  
  cat(cli::bg_cyan(cli::col_red('Done!')))
  cat("\n\n")
  do.call(rbind, Est.trues)
})

library(ggpubr)
out.plots <- list()

plotfn <- function(df){
  do.call(rbind, ET) %>% 
    filter(df == df) %>% 
    tidyr::pivot_longer(`Y.1_.Intercept.`:`Y.1_time`) %>% 
    mutate(
      name = ifelse(grepl("Intercept", name), "b[0]", "b[1]")
    ) %>% 
    tidyr::pivot_wider(names_from = what, values_from = value,
                       names_vary = 'slowest') %>% 
    filter(!(name=='b[0]' & (abs(True) > 10 | abs(Estimate) > 10)),
           !(name=='b[1]' & (abs(True) > 3 | abs(Estimate) > 3))) %>% 
    ggplot() + 
    # geom_abline(intercept = 0, slope = 1, colour = 'black') + 
    # geom_point(aes(x = Estimate, y = True),colour = col, size = .2) + 
    # geom_density(aes(x = Estimate), colour = 'black', trim = TRUE) + 
    # geom_density(aes(x = True), colour = 'red', trim = TRUE) + 
    stat_density(aes(x = Estimate), geom = 'line', colour = 'black', trim = TRUE) + 
    stat_density(aes(x = True), geom = 'line', colour = 'red', trim = TRUE) + 
    labs(x = 'Estimate', y = 'Density', title = bquote(nu==.(df))) + 
    facet_wrap(~name, scales = 'free', labeller = label_parsed) + 
    theme_light() + 
    theme(
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size = 12, face = 'bold'),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
df3 <- plotfn(3)
df5 <- plotfn(5)
df20 <- plotfn(20)
df100 <- plotfn(100)
               
