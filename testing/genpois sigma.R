pgrad <- sapply(1:n, function(i){
  pracma::grad(appxE_GenPoissigma, sigma[[1]], 
               eta = eta[[i]][[1]], Y = dmats$Y[[i]][[1]],
                     tau= tau[[i]][[1]], W = dmats$W[[i]][[1]], w=w, v=v)
})

phess <- lapply(1:n, function(i){
  pracma::hessian(appxE_GenPoissigma, sigma[[1]], 
               eta = eta[[i]][[1]], Y = dmats$Y[[i]][[1]],
               tau= tau[[i]][[1]], W = dmats$W[[i]][[1]], w=w, v=v)
})

mygrad <- sapply(1:n, function(i){
  phi <- dmats$W[[i]][[1]] %*% sigma[[1]]
  Exp <- matrix(0, nrow = length(dmats$Y[[i]][[1]]), ncol = con$gh.nodes)
  for(l in 1:con$gh.nodes){
    Exp[, l] <- w[l] * dmats$Y[[i]][[1]]/(exp(eta[[i]][[1]] + tau[[i]][[1]] * v[l])+dmats$Y[[i]][[1]] * phi)
  }
  Exp <- rowSums(Exp)
  mu <- exp(eta[[i]][[1]])
  a <- (dmats$Y[[i]][[1]] - 1) * Exp - 2*dmats$Y[[i]][[1]]/(1 + phi) + (mu+dmats$Y[[i]][[1]]*phi)/((1+phi)^2)
  crossprod(dmats$W[[i]][[1]] , a)
})

myhess <- lapply(1:n, function(i){
  phi <- dmats$W[[i]][[1]] %*% sigma[[1]]
  Exp <- matrix(0, nrow = length(dmats$Y[[i]][[1]]), ncol = con$gh.nodes)
  for(l in 1:con$gh.nodes){
    Exp[, l] <- w[l] * dmats$Y[[i]][[1]]^2/((exp(eta[[i]][[1]] + tau[[i]][[1]] * v[l])+dmats$Y[[i]][[1]] * phi)^2)
  }
  Exp <- rowSums(Exp)
  mu <- exp(eta[[i]][[1]])
  a <- 3* dmats$Y[[i]][[1]]/((1+phi)^2)-2*(dmats$Y[[i]][[1]]*phi+mu)/((1+phi)^3)-(dmats$Y[[i]][[1]]-1)*Exp
  a <- c(a)
  dmat <- diag(a, nrow = length(dmats$Y[[i]][[1]]))
  t(dmats$W[[i]][[1]]) %*% dmat %*% dmats$W[[i]][[1]]
})

sigma[[1]]-solve(Reduce('+', phess), c(rowSums(pgrad)))
sigma[[1]]-solve(Reduce('+', myhess), c(rowSums(mygrad)))

microbenchmark::microbenchmark(
  `C++ / pracma`= {
    pgrad <- sapply(1:n, function(i){
      pracma::grad(appxE_GenPoissigma, sigma[[1]], 
                   eta = eta[[i]][[1]], Y = dmats$Y[[i]][[1]],
                   tau= tau[[i]][[1]], W = dmats$W[[i]][[1]], w=w, v=v)
    })
    
    phess <- lapply(1:n, function(i){
      pracma::hessian(appxE_GenPoissigma, sigma[[1]], 
                      eta = eta[[i]][[1]], Y = dmats$Y[[i]][[1]],
                      tau= tau[[i]][[1]], W = dmats$W[[i]][[1]], w=w, v=v)
    })
    sigma[[1]]-solve(Reduce('+', phess), c(rowSums(pgrad)))
  },
  `All R, analytical` = {
    mygrad <- sapply(1:n, function(i){
      phi <- dmats$W[[i]][[1]] %*% sigma[[1]]
      Exp <- matrix(0, nrow = length(dmats$Y[[i]][[1]]), ncol = con$gh.nodes)
      for(l in 1:con$gh.nodes){
        Exp[, l] <- w[l] * dmats$Y[[i]][[1]]/(exp(eta[[i]][[1]] + tau[[i]][[1]] * v[l])+dmats$Y[[i]][[1]] * phi)
      }
      Exp <- rowSums(Exp)
      mu <- exp(eta[[i]][[1]])
      a <- (dmats$Y[[i]][[1]] - 1) * Exp - 2*dmats$Y[[i]][[1]]/(1 + phi) + (mu+dmats$Y[[i]][[1]]*phi)/((1+phi)^2)
      crossprod(dmats$W[[i]][[1]] , a)
    })
    
    myhess <- lapply(1:n, function(i){
      phi <- dmats$W[[i]][[1]] %*% sigma[[1]]
      Exp <- matrix(0, nrow = length(dmats$Y[[i]][[1]]), ncol = con$gh.nodes)
      for(l in 1:con$gh.nodes){
        Exp[, l] <- w[l] * dmats$Y[[i]][[1]]^2/((exp(eta[[i]][[1]] + tau[[i]][[1]] * v[l])+dmats$Y[[i]][[1]] * phi)^2)
      }
      Exp <- rowSums(Exp)
      mu <- exp(eta[[i]][[1]])
      a <- 3* dmats$Y[[i]][[1]]/((1+phi)^2)-2*(dmats$Y[[i]][[1]]*phi+mu)/((1+phi)^3)-(dmats$Y[[i]][[1]]-1)*Exp
      a <- c(a)
      dmat <- diag(a, nrow = length(dmats$Y[[i]][[1]]))
      t(dmats$W[[i]][[1]]) %*% dmat %*% dmats$W[[i]][[1]]
    })
    sigma[[1]]-solve(Reduce('+', myhess), c(rowSums(mygrad)))
  }
)
