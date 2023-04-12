ff <- function(df, ...){
  dat <- simData(dof = df, D = diag(c(.25, .05)), beta = t(c(2, 0.33, -0.5, 0.25)),
                 theta = c(-2.9, 0.1),
                 gamma = 0.5, family = list('gaussian'),
                 return.ranefs = T)
  dd <- dat$data
  fit <- joint(list(Y.1 ~ time + cont + bin + (1 + time|id)),
        Surv(survtime, status) ~ bin, data = dd,
        family = list('gaussian'),
        ...)
  return(list(fit = fit,
              REs = dat$ran))
}

test <- ff(2, control = list(verbose = T, tol.rel = 1e-3))
summary(test$fit)



library(cli)

# dfs <- c(2, 3, 4, 5, 6, 7, 8, 9, 10) # change as needed!
dfs <- c(15, 20, 25, 30, 50, 100)
N <- 100

fit.df <- function(df, ...){
  nm <- paste0("Simulation study, df = ", df)
  fn <- paste0('/data/c0061461/fits_tdistn/df_', df, '.RData') # the file name
  cli_progress_bar(name = nm, total = N)
  fits <- true.REs <- vector("list", N)
  for(i in 1:N){
    this <- tryCatch(suppressMessages(ff(df, control = list(tol.rel = 1e-3, return.dmats = FALSE))),
                     error = function(e) NULL)
    if(is.null(this)){
      fits[[i]] <- NULL; true.REs[[i]] <- NULL
    }else{
      fits[[i]] <- this$fit
      true.REs[[i]] <- this$REs
    }
    cli_progress_update()
  }
  cli_progress_done()
  cli_alert_success(sprintf("\nSaving in %s\n", fn))
  out <- list(fits = fits, trueREs = true.REs)
  save(out, file = fn)
}

for(d in dfs) fit.df(d)

