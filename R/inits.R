# ============================================================
# Initial conditions for longitudinal and survival sub-models.
# ============================================================

#' @keywords internal
#' @importFrom pracma nearest_spd
#' @importFrom glmmTMB fixef ranef genpois
Longit.inits <- function(long.formulas, disp.formulas, data, family){
  lapply(long.formulas, function(x) if(!"formula"%in%class(x)) stop('"long.formulas" must be of class "formula"'))
  family.form <- lapply(family, function(f){
    ff <- match.arg(f, c('gaussian', 'binomial', 'poisson', 'genpois', 'Gamma',
                         'negbin'), several.ok = F)
    
    # Set appropriate family ==================
    switch(ff, 
           gaussian = family <- gaussian,
           binomial = family <- binomial,
           poisson = family <- poisson,
           genpois = family <- glmmTMB::genpois(),
           Gamma = family <- Gamma(link='log'),
           negbin = family <- glmmTMB::nbinom2(),
           zip = family <- poisson
    )
    family
  })
  
  K <- length(long.formulas)
  if(K!=length(family)) stop('Uneven long.formulas and family supplied.')
  # Fit using glmmTMB =========================
  fits <- lapply(1:K, function(k){
    fit <- glmmTMB(long.formulas[[k]],
                   family = family.form[[k]], data = data, dispformula = disp.formulas[[k]])
    fit
  })
  
  # Extract ===================================
  # Biomarker names
  markers <- unlist(lapply(1:K, function(k) names(fits[[k]]$modelInfo$respCol)))
  # Fixed-effects
  beta <- do.call(c, lapply(1:K, function(k){
    coef <- glmmTMB::fixef(fits[[k]])$cond
    names(coef) <- paste0(markers[k], '_', names(coef))
    coef
  }))
  # Random effects
  Ds <- lapply(1:K, function(k){
    D <- glmmTMB::VarCorr(fits[[k]])$c$id; dimD <- dim(D)
    D <- matrix(D, nrow = dimD[1], ncol = dimD[2])
    # Checking pos-def. on this D, if not then transform with pracma
    if(any(eigen(D)$values < 0) || (det(D) <= 0)){
      message("Generated covariance matrix not positive semi-definite, occurred for ", markers[k], ".")
      message("This is likely sign of singular fit, continuing anyway...")
      message("\n------- -> Transforming... <- -------\n")
      D <- pracma::nearest_spd(D)
    }
    D
  })
  D <- bDiag(Ds)
  
  # Dispersion ================================
  sigma <- lapply(1:K, function(k){
    if("function"%in%class(family[[k]])) f <- family[[k]]()$family else f <- family[[k]]
    f <- match.arg(f, c('gaussian', 'binomial', 'poisson', 'genpois', 'Gamma', 'genpois', 'negbin'), several.ok = F)
    if(f=='genpois'){ # IDENTITY
      out <- setNames(exp(glmmTMB::fixef(fits[[k]])$disp/2) - 1, paste('phi', k, names(glmmTMB::fixef(fits[[k]])$disp), sep = "_"))
    }else if(f == 'gaussian'){ # VARIANCE
      out <- setNames(glmmTMB::sigma(fits[[k]])^2, paste0('sigma^2_', k))
    }else if(f == 'Gamma'){ # LOG SCALE
      out <- setNames(glmmTMB::fixef(fits[[k]])$disp, paste('shape', k, names(glmmTMB::fixef(fits[[k]])$disp), sep = "_"))
    }else if(f == "negbin"){# LOG SCALE
      out <- setNames(glmmTMB::fixef(fits[[k]])$disp, paste('phi', k, names(glmmTMB::fixef(fits[[k]])$disp), sep = "_"))
    }else{
      out <- 0
    }
    out
  })
  
  # Random effects ============================
  b <- do.call(cbind, lapply(1:K, function(k){
    b <- as.matrix(glmmTMB::ranef(fits[[k]])$cond$id)
    colnames(b) <- paste0(markers[k], '_', colnames(b))
    b
  }))
    
  
  # Populate off-block-diagonal terms in D ====
  off.inds <- which(D == 0, arr.ind = T)
  D[lower.tri(D, F)] <- cov(b)[lower.tri(cov(b), F)]
  D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
    
  list(
    beta.init = beta,
    D.init = D,
    sigma.init = sigma,
    sigma.include = which(unlist(sigma) != 0L),
    b = b,
    responses = markers,
    off.inds = off.inds,
    fits = fits
  )
}

#' @keywords internal
.ToStartStop <- function(data, Tvar){
  this.subj <- list()
  uids <- unique(data$id)
  for(i in uids){
    i.dat <- data[data$id == i, c('id', 'time', Tvar)]
    l <- nrow(i.dat)
    this.one <- cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat[, Tvar]))
    )
    if(this.one[l, 2] != unique(i.dat[, Tvar]))
      rbind(this.one, c(this.one[l, 'id'], this.one[l, 3] + 1e-3, unique(i.dat[, Tvar])))
    this.subj[[i]] <- this.one
  }
  as.data.frame(do.call(rbind, this.subj))
}

#' @keywords internal
.ToRanefForm <- function(time, random.formula, inits.long.frame){
  if(attr(random.formula, 'special') == 'none'){
    out <- model.matrix(as.formula(paste0('~', random.formula)), as.data.frame(time))
  }else if(attr(random.formula, 'special') == 'spline'){
    ww <- which(sapply(lapply(inits.long.frame, class), function(x) "basis"%in%x))
    frame <- inits.long.frame[,ww]
    # out <- model.matrix(as.formula(paste0('~', random.formula)), as.data.frame(time))
    out <- predict(frame, time)
  }
  as.data.frame(out)
}

#' @keywords internal
TimeVarCox <- function(data, b, surv, formulas, b.inds, inits.long){
  # Prepare data
  Tvar <- surv$survtime; Dvar <- surv$status
  ss <- .ToStartStop(data, Tvar); q <- ncol(b) # send to Start-Stop (ss) format
  REs <- as.data.frame(b); REs$id <- 1:nrow(b); K <- length(b.inds)
  ss2 <- merge(ss, REs, 'id')
  invarSurv <- data[,c("id", surv$invar.surv.names, Tvar, Dvar)]
  # invarSurv <- cbind(merge(data[, 'id', drop = F], data.frame(id = 1:surv$n, surv$Smat), 'id'),
  #                    data[,c(Tvar, Dvar)])
  
  ss3 <- merge(ss2, invarSurv, 'id')
  ss3 <- ss3[!duplicated.matrix(ss3), ]
  
  # Create gamma variable
  lhs <- lapply(seq_along(formulas), function(x){
    .ToRanefForm(ss3[,'time1'], formulas[[x]]$random, inits.long$fits[[x]]$frame)
  })
  gamma <- lapply(1:K, function(k){
    unname(rowSums(lhs[[k]] * b[ss3$id, b.inds[[k]]]))
  })
  gamma <- do.call(cbind, gamma)
  colnames(gamma) <- paste0('gamma_', 1:K)
  
  # And join on ...
  ss3 <- cbind(ss3, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3[,Tvar] == ss3$time2, ss3[, Dvar], 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    # paste0('Surv(time1, time2, status2) ~ ', paste0(names(surv$ph$assign), collapse = ' + '), ' + ', paste0('gamma_', 1:K, collapse = ' + '))
    paste0('Surv(time1, time2, status2) ~ ', surv$invar.surv.formula, ' + ', paste0('gamma_', 1:K, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}


#' @keywords internal
parseInits <- function(X, params, inds, inits.long, Omega){
  if(!inherits(X, "list"))
    stop(sQuote("inits"), " must be of class ", sQuote("list"), ".")
  
  poss.names <- c("D", "beta", "sigma", "gamma", "zeta")
  
  check.names <- sapply(names(X), `%in%`, poss.names)
  if(any(!check.names))
    stop(sQuote("inits"), " can only have elements: ", 
         paste0(sapply(poss.names, sQuote), collapse=', '), '.')
  
  # We'll proceed by 'overwriting' the initial conditions
  # obtained via glmmTMB/TimeVarCox. Leave those who aren't specified
  # in the control argument untouched.
  params.new <- params
  Omega.new <- Omega
  if(!is.null(X$D)){
    vD <- vech(X$D)
    wD <- grepl("^D\\[", names(params))
    if(length(vD)!=sum(wD)) # Check dim(D)
      stop("inits$D is improperly dimensioned.")
    params.new[wD] <- vD
    Omega.new$D <- X$D
  }
  if(!is.null(X$beta)){
    lb <- length(unlist(inds$R$beta))
    if(length(X$beta)!=lb)    # Check dim(beta)
      stop("inits$beta is improperly dimensioned.")
    beta.start <- max(grep("^D\\[", names(params))) + 1
    params.new[beta.start:(beta.start + lb - 1)] <- X$beta
    Omega.new$beta <- X$beta
  }
  if(!is.null(X$sigma)){
    if(!inherits(X$sigma, "list"))
      stop("Provide inits$sigma as a list of length K")
    if(length(X$sigma)!=length(inds$R$beta))
      stop("X$sigma must be of length K (set the unused elements to zero).")
    
    il.s <- inits.long$sigma.init
    il.s <- unlist(il.s)
    lS <- length(il.s)
    if(length(unlist(X$sigma))!=lS)
      stop("x$sigma improperly dimensioned.")
    
    to.replace <- names(il.s[il.s != 0])
    new.sigma <- unlist(X$sigma)
    new.sigma <- new.sigma[new.sigma != 0]
    
    # Assume they're in the correct order!
    params.new[match(to.replace, names(params))] <- new.sigma
    Omega.new$sigma <- X$sigma
  }
  if(!is.null(X$gamma)){
    lg <- length(X$gamma)
    if(lg != length(inds$R$beta))
      stop("inits$gamma must be a vector of length K.")
    
    params.new[grepl("^gamma\\_", names(params))] <- X$gamma
    Omega.new$gamma <- X$gamma
  }
  if(!is.null(X$zeta)){
    lz <- length(X$zeta)
    wz <- which(grepl("^zeta\\_", names(params)))
    if(lz != length(wz))
      stop("inits$zeta improperly dimensioned vetor.")
    
    params.new[wz] <- X$zeta
    Omega.new$zeta <- X$zeta
  }
  
  list(params = params.new, Omega = Omega.new)
}


