rm(list=ls())
devtools::load_all('.')
data <- PBC
data <- na.omit(data[,c('id', 'time', 'drug', 'platelets', 'albumin','survtime', 'status')])
data$temp <- c(scale(data$platelets))

long.formulas1 <- list(
  temp ~ time * drug + (1 + time|id)
)
long.formulas2 <- list(
  platelets ~ time * drug + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ drug

gauss.fit <- joint(long.formulas1,
                   surv.formula,
                   data = data,
                   family = list('gaussian'))

poiss.fit <- joint(long.formulas2,
                   surv.formula,
                   data = data,
                   family = list('poisson'), control = list(verbose = T))

gp.fit <- joint(long.formulas2,
                surv.formula,
                data = data,
                family = list('genpois'))

ROCs <- setNames(vector('list', 3), paste0('Tstart = ', c(3,6,8)))
for(TT in c(3,6,8)){
  out <- setNames(vector('list', 3), c('Gaussian', 'Poisson', 'genPois'))
  for(f in c('Gaussian', 'Poisson', 'genPois')){
    if(f == 'Gaussian') 
      out[[f]] <- ROC(gauss.fit, data, TT, 2, control = list(nsim = 25))
    else if(f == 'Poisson')
      out[[f]] <- ROC(poiss.fit, data, TT, 2, control = list(nsim = 25))
    else 
      out[[f]] <- ROC(gp.fit, data, TT, 2, control = list(nsim = 25))
    
    message('TT: ', TT, ' done for family ', f)
    cat('\n')
  }
  ROCs[[paste0('Tstart = ', TT)]] <- out
}


# (3, 5]
gauss.roc <- ROC(gauss.fit, data, 8, 2, control = list(nsim = 25))
poiss.roc <- ROC(poiss.fit, data, 8, 2, control = list(nsim = 25))
gp.roc <-    ROC(gp.fit,    data, 8, 2, control = list(nsim = 25))

# (6, 8]

# (8, 10]


# Unrelated: Ensuring Omega.draw works ------------------------------------

test.fit <- joint(long.formulas = list(albumin ~ time + (1+  time|id),  platelets ~ time * drug + (1 + time|id)), 
                  surv.formula = surv.formula, data = data, family = list('gaussian', 'poisson'))
