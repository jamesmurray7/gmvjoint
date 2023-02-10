load('~/Downloads/MVNsimsets.RData')

ff <- function(x, f){
  if(f == 'binomial')
    long.formulas <- list(Y.1 ~ time + cont + bin + (1|id))
  else
    long.formulas <- list(Y.1 ~ time + cont + bin + (1 + time|id))
  
  out <- tryCatch(suppressMessages(joint(
    long.formulas, 
    Surv(survtime, status) ~ bin,
    data = x, # x single data set!
    family = list(f)
  )),
  error = function(e) NULL)
  
  out
}

nm <- names(simsets)
fitsets <- setNames(vector("list", length(nm)),
                    nm)
for(s in seq_along(nm)){
  snm <- nm[s]
  fam <- trimws(stringr::str_extract(snm, "\\w+$"))
  thisxlist <- simsets[[s]]
  out <- vector('list', 100)
  pb <- utils::txtProgressBar(max=100, style=3)
  for(j in 1:100){
    out[[j]] <- ff(thisxlist[[j]], fam)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)
  fitsets[[s]] <- out
  cat(sprintf("\n%d out of %d done.\n", s, length(nm)))
}