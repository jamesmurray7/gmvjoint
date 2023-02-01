# Updated version of parseCoxph
parseCoxph2 <- function(surv.formula, data, center = TRUE){
  survdata <- data[!duplicated(data[, 'id']), ]; n <- nrow(survdata)
  ph <- coxph(surv.formula, survdata, x = T)
  
  survdata <- data.frame(id = survdata$id, ph$x, survtime = ph$y[, 1], status = ph$y[, 2])
  # One row per id for each survival covariate.
  Smat <- model.matrix(ph)
  pS <- ncol(Smat)
  # Use survfit.coxph with newdata argument to avoid warnings.
  sf <- survfit(ph, newdata = as.data.frame(cbind(Smat, survdata[, c("survtime", "status")])))
  ft <- sf$time[sf$n.event >= 1]     # failure times
  nev <- sf$n.event[sf$n.event >= 1] # Number of failures per failure time
  
  # Survival indicators
  Delta <- as.list(survdata$status)
  
  # Scale x variables
  if(center){ 
    survdata[,2:(pS+1)] <- apply(survdata[, 2:(pS+1), drop = F], 2, scale, scale=F)
    Smat <- apply(Smat, 2, scale, scale = F)
  }
  
  # Names used in surv.formula (for use in other functions).
  parsed.ph <- extract.surv.process(ph)
  survtime <- parsed.ph$Time
  status <- parsed.ph$Status
  
  # Return ----
  structure(list(
    survdata = survdata, ph = ph, 
    Smat = Smat, Delta = Delta,
    n = n, ft = ft, nev = nev,
    survtime = survtime, status = status
  ), class = 'parseCoxph')
  
}

print.parseCoxph <- function(x, ...){
  if(!inherits(x, 'parseCoxph')) stop("Only for use with objects of class 'parseCoxph'.\n")
  print(x$ph)
  invisible(x)
}

surv.mod2 <- function(surv, formulas, l0.init){
  # unpack parseCoxph object
  n <- surv$n; K <- length(formulas); l0 <- l0.init
  ft <- surv$ft; survtime.name <- surv$survtime; status.name <- surv$status
  
  # Time-invariant covariates
  S <- lapply(1:n, function(i) as.matrix(surv$Smat[i,,drop=F]))
  
  # vector function of time to use in hazard
  Wk <- lapply(1:K, function(k) formulas[[k]]$random)
  
  # Fu, design matrix of _all_ failure times.
  ft.df <- data.frame(time = ft)
  Fu.all <- do.call(cbind, lapply(Wk, function(f) model.matrix(as.formula(paste0("~", f)), ft.df)))
  
  # Failure times and status list for each id = i,...,n.
  TiDi <- lapply(1:n, function(x) surv$survdata[surv$survdata$id == x, c('survtime', 'status')])
  
  # Create Fi,
  Fi <- lapply(1:n, function(i){
    T.df <- data.frame(time = TiDi[[i]]$survtime)
    do.call(cbind, lapply(Wk, function(f) model.matrix(as.formula(paste0("~", f)), T.df)))
  })
  
  # Populate other items
  Fu <- l0u <- l0i <- surv.times <- SS <- vector('list', n)
  for(i in 1:n){
    Ti <- TiDi[[i]]$survtime; Di <- TiDi[[i]]$status
    # Failure times survived (up-to-and-including Ti).
    surv.times[[i]] <- which(ft <= Ti) # Store indiced
    St <- ft[surv.times[[i]]]          # The actual times
    if(length(St)){   # Design matrices of 
      Fu[[i]] <- Fu.all[surv.times[[i]], , drop = F]
      l0u[[i]] <- l0[surv.times[[i]]]
    }else{            # Case when individual is censored before first failure time.
      Fu[[i]] <- matrix(0, nrow = 1, ncol = ncol(Fu.all))
      l0u[[i]] <- 0
    } 
    SS[[i]] <- apply(S[[i]],2,rep,nrow(Fu[[i]]))
    if(!"matrix"%in%class(SS[[i]])) SS[[i]] <- t(SS[[i]])
    if(Di == 1L) l0i[[i]] <- l0[which(ft == Ti)] else l0i[[i]] <- 0
  }
  
  # Return
  list(
    ft = ft, ft.mat = Fu.all, nev = surv$nev, surv.times = surv.times,
    l0 = l0, l0i = l0i, l0u = l0u, 
    Fi = Fi, Fu = Fu, Tis = sapply(TiDi, function(x) c(x[1]$survtime), simplify = T),
    S = S, SS = SS, q = ncol(Fu.all)
  )
  
}
