devtools::load_all('.')
df <- survival::pbcseq[1:10,]
library(Rcpp)
library(RcppArmadillo)
sourceCpp('dev/cpp-modelmatrices.cpp')

CreateFixedMatrices(ff, df, unique(df$id))


d <- simData()$data

df <- d
formulae <- list(
  ~ bin * time + I(time^2) + cont,
  ~ time + cont + bin
)

X <- CreateFixedMatrices(formulae, df, 1:250)
