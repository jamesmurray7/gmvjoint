cat("\f")
message("Reloading---\nDevelopment version---\n")
sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}

sourceC <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "cpp$")) {
    if(trace) cat(nm,":")
    sourceCpp(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}
library(glmmTMB)
library(survival)
library(Rcpp)
library(RcppArmadillo)
rm_ <- function() rm(list = setdiff(ls(), c("sourceDir", "sourceC", "rm_")))
sourceDir("dev/DEVELOPMENT/R/")
sourceC("dev/DEVELOPMENT/src/")
cat("ready!\n\n")
