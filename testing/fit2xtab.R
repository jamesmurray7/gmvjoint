fit2xtab <- function(fit, max.row = NULL){
  if(is.null(fit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3)
  
  # Model fit info
  K <- length(fit$ModelInfo$ResponseInfo)
  responses <- lapply(sapply(fit$ModelInfo$ResponseInfo, strsplit, '\\s\\('), el, 1)
  families <- unlist(fit$ModelInfo$family)
  # Standard errors and parameter estimates.
  SE <- fit$SE
  D <- fit$co$D
  betas <- fit$co$beta
  sigmas <- unlist(fit$co$sigma)
  gammas <- fit$co$gamma
  zetas <- setNames(fit$co$zeta, paste0('zeta_', 1:length(fit$co$zeta)))
  
  MakeTables <- lapply(1:K, function(k){
    
    nb <- names(betas)[grepl(responses[[k]], names(betas))]
    nb2 <- paste0('beta_{', k, (seq(0, (length(nb) - 1))), '}')
    beta <- setNames(betas[grepl(responses[[k]], names(betas))], nb2)
    sigma <- setNames(sigmas[k], 
                      if(families[k] == "genpois")
                        paste0("phi_", k)
                      else if(families[k] == 'Gamma')
                        paste0("shape_", k)
                      else if(families[k] == 'gaussian')
                        paste0('sigma^2_', k)
                      else
                        'NO DISP')
    
    betak <- c(beta)
    if(sigma != "NO DISP") betak <- c(betak, sigma)
    parameter <- names(betak)
    
    rSE <- SE[match(names(betak), names(SE))]#SE associated with these coeffs
    
    gamma <- setNames(gammas[k], paste0('gamma_', k))
    
    kk <- c(beta, sigma, gamma)
    kk.names.lookup <- c(nb, names(sigma), names(gamma))
    
    kSE <- SE[match(kk.names.lookup, names(SE))]#SE associated with these coeffs
    
    lb <- kk - qz * kSE; ub <- kk + qz * kSE
    
    this.out <- setNames(data.frame(.to3dp(kk), .to3dp(kSE), .to3dp(lb), .to3dp(ub)),
                         c('Estimate', 'SE', '2.5%', '97.5%'))
    this.out
  })
  
  tab <- do.call(rbind, MakeTables)
  # Append zeta terms to bottom, time invariant so report separately
  SEz <- SE[grepl('^zeta', names(SE))]; lb <- zetas - qz * SEz; ub <- zetas + qz * SEz
  zet <- data.frame(.to3dp(zetas), .to3dp(SEz), .to3dp(lb), .to3dp(ub))
  names(zet) <- names(tab)
  tab <- rbind(tab, zet)
  nr <- nrow(tab)
  
  tab$Parameter <- rownames(tab)
  tab <- tab[!grepl('NO DISP', tab$Parameter),]
  tab$Parameter <- paste0('$\\', tab$Parameter, '$')
  
  tab2 <- as.data.frame(cbind(Parameter = tab$Parameter, apply(tab[, -5], 2, function(x) format(round(x, 3), nsmall = 3))))
  tab3 <- cbind(Parameter = tab$Parameter, `Mean (SE)` = paste0(tab2$Estimate, ' (', tab2$SE, ')'), 
                `95% CI` = paste0('[', tab2$`2.5%`, ', ', tab2$`97.5%`, ']'))
  
  # Splitting out into multiple columns -->
  if(nr > 15 && is.null(max.row)){
    cat('Consider breaking at a certain number of rows and presenting a "wider" table.\nnrows: ', nr, '\n') 
  }
  
  if(!is.null(max.row)){
    if(nr <= max.row) stop('max.row must exceed the number of rows in output table: ', nr, '.\n')
    
    # Work out how many 'cbinds' we'll need to do.
    num.splits <- nr %/% max.row
    nr.each <- suppressWarnings(split(1:nr, rep(1:num.splits, each = ceiling(nr/num.splits))))
    split.tab <- lapply(nr.each, function(k){
      x <- tab3[k,]
      while(nrow(x) < max.row){
        x <- rbind(x, c('-','-','-'))
      }
      x
    })
    
    tab3 <- do.call(cbind, split.tab)
  }
  
  xt <- xtable::xtable(tab3,
                       caption = paste0("Elapsed time for approximate EM algorithm to converge and SE calculation was ", round(fit$elapsed.time['EM time'] + 
                                                                                                                               fit$elapsed.time["Post processing"], 2), " seconds."))
  
  print(xt,
        include.rownames = FALSE,
        sanitize.text.function = identity)
  
}

