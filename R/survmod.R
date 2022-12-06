# Survival objects


#' Parsing the survival formula and constructing all survival-related data objects.
#' @param surv.formula A formula readable by `coxph`.#'
#' @param data a set of data containing covariate information for variables
#'   named by `surv.formula`. Can be of any 'completeness', as the function 
#'   returns a reduced set.
#'
#' @returns A list containing 
#' 
#' * `survdata`: Reduced version of `data`, with only one row per subject, with covariates 
#'  specified by `surv.formula` along with survival time and failure status. 
#' * `ph`: model fit from `coxph`. 
#' * `n`: Number of unique subjects. 
#' * `Delta`: List of failure indicators for each subject (1=failed). 
#' 
#' @export
#' @examples 
#' 
#' data = simData()$data
#' parseCoxph(Surv(survtime, status) ~ bin, data = data)
parseCoxph <- function(surv.formula, data){
  survdata <- data[!duplicated(data[, 'id']), ]; n <- nrow(survdata)
  ph <- coxph(surv.formula, survdata, x = T)
  
  survdata <- data.frame(id = survdata$id, ph$x, survtime = ph$y[, 1], status = ph$y[, 2])
  pS <- ncol(ph$x)
  sf <- survfit(ph)
  ft <- sf$time[sf$n.event >= 1]     # failure times
  nev <- sf$n.event[sf$n.event >= 1] # Number of failures per failure time
  
  Delta <- as.list(survdata$status)
  
  # Return ----
  list(
    survdata = survdata, ph = ph, n = n, Delta = Delta
  )
}


# Create survival data objects based on random effects formula(e), a ph fit,
# the survival data and an initial estimation of \eqn{\lambda_0}.
#' @keywords internal
surv.mod <- function(ph, survdata, formulas, l0.init){
  uids <- unique(survdata$id); n <- length(uids); K <- length(formulas)
  l0 <- l0.init; splines <- F
  
  # initialise empty stores
  l0i <- vector('numeric', n)
  l0u <- surv.times <- vector('list', n)
  
  # Failure times
  ft <- coxph.detail(ph)$time
  
  # First loop over subjects to create K-invariant objects.
  for(i in as.numeric(uids)){
    # slice this id
    survtime <- survdata[survdata$id == i, 'survtime']
    status <- survdata[survdata$id == i, 'status']
    # INDICES of survived times
    surv.times[[i]] <- which(ft <= survtime)
    # Fu, and l0u
    st <- ft[which(ft <= survtime)] # survived times
    if(length(st) > 0){
      l0u[[i]] <- l0[which(ft <= survtime)]
    }else{ # Those who are censored before first failure time
      l0u[[i]] <- 0
    }
    if(status == 1) l0i[i] <- l0[which(ft == survtime)] else l0i[i] <- 0
  }
  
  # Second loop over formulas and subjects to create i x K object. ----
  
  sv <- lapply(formulas, function(formulas){
    if((!is.null(attr(formulas$random, 'special'))) & attr(formulas$random, 'special') == 'spline'){
      splines <- T
      survdata$time <- unname(ph$y[,1])
      newSurvdata <- as.data.frame(cbind(id = survdata$id, model.matrix(as.formula(paste0('~', formulas$random)), survdata)))
      q <- ncol(newSurvdata) - 1
      .getFi <- function(time, q = q){
        id <- which(survdata$time == time)[1] # Take the first one in case of ties -- same design matrix regardless
        as.matrix(newSurvdata[id, -1, drop = F])
      } 
      .getFu <- function(times, q = q){
        as.matrix(newSurvdata[match(times, survdata$time), -1, drop = F])
      }
      Fi <- lapply(survdata$survtime, .getFi)
      ftmat <- .getFu(ft)
    }else{ # Case 2:: Anything else // Maybe more specials in future(?)
      random <- formulas$random
      q <- length(el(strsplit(random, '\\+|\\*|\\:')))
      # Define two functions to construct F_i and F_u ----
      .getFi <- function(time, q = q){
        Fi <- matrix(NA, nrow = 1, ncol = q)
        for(i in 1:q) Fi[, i] <- time^(i - 1)
        Fi
      }
      .getFu <- function(times, q = q){
        out <- sapply(1:q, function(i) times ^ (i - 1))
        if(!"matrix"%in%class(out)) out <- t(out)
        out
      }
      # Generate F_i and ftmat
      Fi <- lapply(survdata$survtime, .getFi, q)
      ftmat <- .getFu(ft, q)
    }
    
    # loop over subjects -----
    Fu <- vector('list', n)
    for(i in as.numeric(uids)){
      # slice this id
      survtime <- survdata[survdata$id == i, 'survtime']
      status <- survdata[survdata$id == i, 'status']
      # INDICES of survived times
      surv.times[[i]] <- which(ft <= survtime)
      # Fu, and l0u
      st <- ft[which(ft <= survtime)] # survived times
      if(length(st) > 0){
        Fu[[i]] <- .getFu(st, q)
      }else{ # Those who are censored before first failure time
        Fu[[i]] <- do.call(cbind, replicate(q, 0, simplify = F))
      }
    }
    list(Fu = Fu, Fi = Fi, ftmat = ftmat)
  })
  
  # Collate Fi, Fu, ftmat.
  Fi <- lapply(lapply(1:n, function(i){
    lapply(1:K, function(k){
      sv[[k]]$Fi[[i]]
    })
  }), function(x) do.call(cbind, x))
  
  Fu <- lapply(lapply(1:n, function(i){
    lapply(1:K, function(k){
      sv[[k]]$Fu[[i]]
    })
  }), function(x) do.call(cbind, x))
  
  ftmat <- lapply(1:K, function(k){
    sv[[k]]$ftmat
  })
  
  # Design matrix SS and rowvec S
  S <- lapply(1:n, function(i) ph$x[i, , drop = F])
  SS <- lapply(1:n, function(i){
    out <- apply(S[[i]], 2, rep, nrow(Fu[[i]]))
    if("numeric"%in%class(out)) out <- t(out)
    out
  })
  
  # Return list ----
  return(list(
    ft = ft, ft.mat = ftmat, nev = coxph.detail(ph)$nevent, surv.times = surv.times,
    l0 = l0, l0i = as.list(l0i), l0u = l0u, Fi = Fi, Fu = Fu, Tis = survdata$survtime,
    S = S, SS = SS, q = q
  ))
}