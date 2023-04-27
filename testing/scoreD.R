D.update.bh <- Map(function(S, b){
  SqrtDiagS <- sqrt(diag(S))
  SqrtS <- t(chol(S))
  store <- matrix(0, ncol = q, nrow = q)
  for(j in 1:q) store[j,j] <- sum(w * ((b[j] + SqrtDiagS[j]*v)^2))
  for(p in 1:(q-1))
    for(r in (p+1):q)
      for(l1 in 1:length(w))
        for(l2 in 1:length(v))
          store[p,r] <- w[l1]*w[l2]*(b[p]+SqrtS[p,c(p,r)]%*%c(v[l1], v[l2])) * 
                          (b[r]+SqrtS[r,c(p,r)]%*%c(v[l1], v[l2])) + store[p,r]
  0.5 * (store + t(store))
}, S = Sigma, b = b.hat)

pmult <- diag(rep(1, q))
pmult[pmult==0] <- .5
D.update2 <- mapply(function(Sigma, b){
  out <- Sigma + tcrossprod(b)
  out * pmult
} , Sigma = Sigma, b = b.hat, SIMPLIFY = F)

Reduce('+', D.update)/n #-
Reduce('+', D.update2)/n #-
Reduce('+', D.update.bh)/n

ScoreD <- function(D, S, b){
  Dinv <- solve(D)
  vech(t(.5 * Dinv %*% (S + tcrossprod(b))  %*% Dinv - .5 * Dinv))
}

InfD <- function(D, S, b){
  vD <- vech(D)
  XX <- pmax(abs(vD), 1)
  out <- matrix(0, length(vD), length(vD))
  S0 <- ScoreD(D, S, b)
  for(d in 1:length(vD)){
    vD.d <- vD
    vD.d[d] <- vD[d] + 1e-3 * XX[d]
    thisD <- vech2mat(vD.d, q)
    diffS <- ScoreD(thisD, S, b) - S0
    out[,d] <- diffS / (vD.d[d]-vD[d])
  }
  out
}

ScoreDs <- mapply(function(b, S) ScoreD(D, S, b), S = Sigma, b = b.hat)
InfDs <- Map(function(b, S) InfD(D, S, b), S = Sigma, b = b.hat)

vech2mat(vech(D) - solve(Reduce('+', InfDs), rowSums(ScoreDs)), q) # mathces most closesly



