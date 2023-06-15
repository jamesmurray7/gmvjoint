#' Print an LaTeX-ready \code{xtable} for a \code{joint} object.
#' 
#' @description Prints an \code{xtable} output for a fitted \code{joint} object to the console,
#' or to a specified save location
#'
#' @param x a joint model fit by the \code{joint} function.
#' @param caption character, specifies the \code{caption} argument of \code{xtable}. By default
#' this takes value \code{NULL}, which results in a generic caption being generated.
#' @param label character, specifies the \code{label} argument of \code{xtable}.
#' @param align character, specifies the \code{align} argument of \code{xtable}. Note by default
#' this is \code{NULL}, as alignment is done internally.
#' @param digits integer, specifies the \code{digits} argument of \code{xtable}. Note by default
#' this is \code{NULL}, as argument \code{dp} controls this (but can be specified through this, 
#' too).
#' @param display character, specifies the \code{display} argument of \code{xtable}.
#' @param auto logical, specifies the \code{auto} argument of \code{xtable}. Defaults to 
#' \code{FALSE}. Not recommended to change.
#' @param p.val logical, should p-values be returned? Defaults to \code{p.val = FALSE}.
#' @param max.row integer, the number of rows after which the table is `broken' vertically
#' and merged horizontally; useful for long tables. Defaults to \code{max.row = NULL} which
#' results in one long table. Note that this can be quite finicky, so trial and error may be 
#' required.
#' @param dp integer, the number of decimal places to round the estimate, standard error and
#' confidence intervals to; defaults to \code{dp = 3}.
#' @param vcov logical, should the half-vecorisation of the block diagonal of covariance 
#' matrix be reported? Default is \code{vcov = FALSE}.
#' @param capture logical, should the printed \code{xtable} output be saved anywhere instead
#' of just printed to the console? Defaults to \code{capture = FALSE}.
#' @param capture.location character, if \code{capture = TRUE}, this should specify what
#' \emph{file} it should be saved in. Defaults to \code{capture.location = ""}.
#' @param hlines character, specifies which horizontal lines are used in the outputted
#' LaTeX table. Supply a character string which contains \code{"top"}, \code{"middle"} and/or
#' \code{"bottom"} (in any order) to specify a \code{toprule}; \code{midrule} and 
#' \code{bottomrule} in the table. If \code{booktabs = FALSE}, then these will simply 
#' be \code{hline}s. For instance \code{hlines = "top-middle-bottom"} prints all three;
#' whilst \code{hlines = "middle-bottom"} skips the \code{toprule}.
#' @param booktabs logical, if \code{booktabs = TRUE} (the default) then \code{toprule};
#' \code{midrule} and \code{bottomrule} replace the usual \code{hline}s.
#' @param size character, LaTeX size to be placed before the tabular environment, defaults
#' to \code{size = "footnotesize"}; replace with \code{"normalsize"} if wanted.
#' @param ... additional arguments, none used.
#'
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @method xtable joint
#' @return A LaTeX-ready \code{xtable} print-out of the joint model. A list containing 
#' constituent tables is also returned invisibly, along with the final \code{xtable} output.
#' @seealso \code{\link{joint}}
#' @export
#'
#' @examples
#' # Bivariate joint model ------------------------------------------------
#' require(xtable)
#' data <- simData(n = 100)$data
#' long.formula <- list(
#'   Y.1 ~ time + cont + bin + (1 + time|id),
#'   Y.2 ~ time + cont + bin + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ cont + bin
#' family <- list("gaussian", "gaussian")
#' fit <- joint(long.formula, surv.formula, data, family)
#' xtable(fit)
#' # Example of arguments: add dummy caption, add p-values.
#' xtable(fit, p.val = TRUE, dp = 4, caption = "This is a caption")
#' # Change size, place horizontal lines everywhere
#' xtable(fit, size = "normalsize", hlines = c("top-middle-bottom"))
#' # Make a wider table without booktabs 
#' xtable(fit, booktabs = FALSE, max.row = 6)
xtable.joint <- function(x, caption = NULL, label = NULL, align = NULL, digits = NULL, # Default arguments
                         display = NULL, auto = FALSE,                                 # (Needed for S3 consistency)
                         p.val = FALSE, max.row = NULL, dp = 3,
                         vcov = FALSE, capture = FALSE, capture.location = "", 
                         hlines = "middle-bottom", booktabs = TRUE, size = "footnotesize",
                         ...){
  # Checks
  if(!inherits(x, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(x$SE)) stop('re-run with post.process = TRUE')
  if (!requireNamespace("xtable", quietly = TRUE)) 
    stop("'xtable' is required.\n")
  # Random stuff
  qz <- qnorm(.975)
  if(!is.null(digits) & is.null(dp)) dp <- digits # Safety in case digits supplied for some reason.
  .toXdp <- function(x) format(round(x, dp), nsmall = dp)
  
  # Model fit info
  K <- x$ModelInfo$K
  responses <-x$ModelInfo$Resps
  families <- unlist(x$ModelInfo$family)
  # This will be based on summary
  s <- summary(x)
  L <- s$Longits; S <- s$Survs
  # Put zeta(s) at end
  S.gamma <- S[(1+length(x$coeffs$zeta)):nrow(S),]
  S.zeta <- S[1:length(x$coeffs$zeta),]
  
  # If (block diagonal) of vcov elements are required, form these now.
  if(vcov){
    D <- x$coeffs$D
    SE.vD <- x$SE[grepl("^D\\[", names(x$SE))]
    vD <- setNames(vech(D), names(SE.vD))
    
    b.inds <- x$ModelInfo$inds$R$b
    Dtabs <- lapply(seq_along(b.inds), function(i){
      x <- b.inds[[i]]
      Dx <- D[x,x]; vDx <- vech(Dx)
      inds <- which(vD %in% vDx)
      xSE <- SE.vD[inds]
      
      Parameter <- paste0("D_{", i,",",apply(which(lower.tri(Dx, T), arr.ind = T) - 1, 1, paste, collapse=''),"}")
      MSE <- paste0(.toXdp(vDx), " (", .toXdp(xSE), ")")
      lb <- vDx - qz * xSE; ub <- vDx + qz * xSE
      CI <- paste0('[', .toXdp(lb),', ',.toXdp(ub),']')
      
      tab.Dx <- cbind(Parameter, MSE, CI)
      if(p.val)
        return(cbind(tab.Dx, rep('{}', length(vDx))))
      else
        return(tab.Dx)
    })
  }
  
  # Rearrange to get {longit, gamma} for each k=1,dots,K response
  RespChunks <- lapply(1:K, function(k){
    Lk <- L[[k]];
    Respk <- responses[[k]]
    Sk <- S.gamma[k,]
    df <- rbind(Lk, Sk)
    # Get parameter names
    Parameter <- row.names(df)
    # Rename gamma to gamma_k
    Parameter <- ifelse(grepl("gamma", Parameter), gsub(paste0("\\_", Respk), paste0("_",k), Parameter),
                        Parameter)
    # Rename beta to beta_{k(0:Pk)}
    fixefs.rhs <- which(grepl(Respk, Parameter)) - 1
    Parameter <- ifelse(grepl(Respk, Parameter), 
                         paste0("beta_{", k, fixefs.rhs, "}"),
                        Parameter)
    MSE <- paste0(.toXdp(df$Estimate), " (", .toXdp(df$SE), ")")
    CI <- paste0('[', .toXdp(df$`2.5%`),', ',.toXdp(df$`97.5%`),']')
    out <- cbind(Parameter, MSE, CI)
    if(p.val){
      p <- ifelse(df$`p-value` < 1e-3, "< 0.001", 
                  format(round(df$`p-value`, 3), nsmall = 3, justify = 'right', width = 7))
      out <- cbind(out, p)
    }
    if(vcov) out <- rbind(Dtabs[[k]], out)
    out
  })
  
  tab <- do.call(rbind, RespChunks)
  # zeta is time invariant so report separately at foot of table
  zeta.names <- row.names(S.zeta)
  if(length(zeta.names) > 1) Parameter <- paste0("zeta_", 1:length(zeta.names)) else Parameter <- 'zeta'
  MSE <- paste0(.toXdp(S.zeta$Estimate), " (", .toXdp(S.zeta$SE), ")")
  CI <- paste0('[', .toXdp(S.zeta$`2.5%`),', ',.toXdp(S.zeta$`97.5%`),']')
  tab.zeta <- cbind(Parameter, MSE, CI)
  if(p.val){
    p <- ifelse(S.zeta$`p-value` < 1e-3, "< 0.001", 
                format(round(S.zeta$`p-value`, 3), nsmall = 3, justify = 'right', width = 7))
    tab.zeta <- cbind(tab.zeta, p)
  }

  tab2 <- as.data.frame(rbind(tab, tab.zeta), stringsAsFactors = FALSE)
  tab2$Parameter <- paste0('$\\', tab2$Parameter, '$')
  
  # Splitting out into multiple columns -->
  nr <- nrow(tab2)
  if(nr > 15 && is.null(max.row)){
    cat('Consider breaking at a certain number of rows and presenting a "wider" table.\nnrows: ', nr, '\n') 
  }
  names(tab2)[2] <- "Estimate (SE)"; names(tab2)[3] <- "95\\% CI"
  if(p.val) names(tab2)[4] <- "p-value"
  tab3 <- tab2
  align <- if(!p.val) "cl|rr" else "cl|rrr"
  if(!is.null(max.row)){
    # How many rows should be in each split?
    num.splits <- ceiling(nr / max.row)
    split.indices <- suppressWarnings(split(1:nr, 
                           rep(1:num.splits, each = max.row)))
    split.indices <- lapply(split.indices, function(x){
      x[which(x<=nr)]
    })
    split.tabs <- lapply(split.indices, function(k){
      this <- tab3[k,]
      .nr <- nrow(this)
      while(.nr < max.row){
        this <- rbind(this, rep('{}', ncol(tab3)))
        .nr <- nrow(this)
      }
      this
    })
    tab3 <- do.call(cbind, split.tabs)
    names(tab3) <- rep(names(tab2), num.splits)
    rep.align <- paste0(rep(gsub("^c", "", align), num.splits-1), collapse = '')
    align <- paste0(align, rep.align, collapse = '')
    align <- gsub("\\|", "", align)
  }

  if(is.null(caption)){
    resps <- unlist(responses)
    if(K >= 2){
      resps2 <- paste0(resps[-K], collapse=', ')
      report.resps <- paste0(resps2, ', and ', resps[K])
    }else if(K == 2){
      report.resps <- paste0(resps[1], " and ", resps[2])
    }else{
      report.resps <- resps[1]
    }
    
    caption <- sprintf("Parameter estimates (SE: standard error) for the joint analysis of %s. Elapsed time for the approximate EM algorithm to converge and SE calculation was %.3f seconds. Total Computation time was %.3f seconds.",
                       report.resps, x$elapsed.time['EM time'] + x$elapsed.time["Post processing"], x$elapsed.time["Total Computation time"])
  }else{
    caption <- caption
  }
  # Work out hline input
  hline.after <- c(-1, 0, nrow(tab3))
  if(!grepl("top", hlines)) hline.after[1] <- NA
  if(!grepl("middle", hlines)) hline.after[2] <- NA
  if(!grepl("bottom", hlines)) hline.after[3] <- NA
  hline.after <- hline.after[!is.na(hline.after)]
  
  # Create the xtable
  xt <- xtable::xtable(tab3, caption = caption, label = label, align = align, 
                       digits = digits, display = display, auto = auto)
  
  if(capture){
    if(nchar(capture.location)==0L) stop("Please specify capture.location.")
    capture.output(print(xt,
                         include.rownames = FALSE,
                         sanitize.text.function = identity,
                         booktabs = booktabs,
                         size = size,
                         hline.after = hline.after),
                   file = capture.location)
    cat(sprintf("Saved in %s\n", capture.location))
  }else{
    print(xt, 
          include.rownames = FALSE,
          sanitize.text.function = identity,
          booktabs = booktabs,
          size = size,
          hline.after = hline.after)
  }
  
  return(invisible(list(
    Dtabs = if(vcov) Dtabs else NULL,
    RespChunks = RespChunks,
    zeta = tab.zeta,
    xtab = xt
  )))
}

