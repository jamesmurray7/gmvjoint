eta
t
t2 <- .5 * t^2
t2

true <- exp(eta)

exp(eta+t2)

gherm <- matrix(0,10,3)
for(l in 1:3) gherm[,l] <- w[l] * exp(eta + v[l] * t)
rowSums(gherm)



# '''
true <- log(1 + exp(eta))
appx <- log(1+exp(eta + t2)) - (exp(t2) - 1) * exp(eta + t2)/(2*(1+exp(eta + t2))^2)
true
appx


gherm <- matrix(0,10,3)
for(l in 1:3) gherm[,l] <- w[l] * log(exp(eta + v[l] * t)+1)
rowSums(gherm)




