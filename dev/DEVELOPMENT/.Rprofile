cat("\f")  
message("Development version of gmvjoint")
options(prompt = crayon::green("Dev> "))
library(Rcpp)
library(RcppArmadillo)
suppressWarnings(suppressMessages(library(glmmTMB)))
library(survival)

scc <- function(){
  cat("\nfuns.cpp")
  sourceCpp('src/funs.cpp')
  cat("\nGP1.cpp")
  sourceCpp("src/GP1_PMF_SCALAR.cpp")
  cat("\ndynpreds.cpp")
  sourceCpp("src/dynpreds.cpp")
  cat("\n")
}
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

scR <- function() sourceDir('R/')
message("\n Loading R/ directory\n")
scR()


message("Loading src/ directory\n")
scc()

message("Ready!\n")
