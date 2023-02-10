all.sims <- expand.grid(mi = c(5, 10, 15), 
                        n = 250,
                        family = c('gaussian', 'poisson', 'binomial'))

D <- diag(c(0.25, 0.05))

# Fixed effects
beta <- t(c(2, -0.1, 0.1, 0.2))

# gamma & zeta
zeta <- c(-0.0, 0.3)

.sim <- function(n, mi, f){
  if(f == 'binomial'){
    D <- matrix(0.25, 1, 1) 
    random.formula <- list(~1)
  }else{
    D <- diag(c(0.25, 0.05))
    random.formula <- NULL
  }
  if(f != 'poisson') gamma <- c(-0.50) else gamma <- c(0.50)
  
  simData(n = n, ntms = mi, family = as.list(f),
          beta = beta, sigma = c(0.16), D = D,
          gamma = gamma, zeta = zeta, 
          random.formula = random.formula,
          theta = c(-2.5, 0.1))$data
}

simsets <- setNames(apply(all.sims, 1, 
                 function(x) replicate(100, .sim(as.numeric(x[2]), 
                                                 as.numeric(x[1]), x[3]),
                                       simplify = F)),
                 apply(all.sims, 1, function(x) 
                   paste0('n = 250, mi = ', as.numeric(x[1]), ', fam = ', (x[3]))))

save(simsets, file = '~/Downloads/MVNsimsets.RData')
