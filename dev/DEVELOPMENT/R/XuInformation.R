ScoreOmega <- function(Omega.perturbed, dmats, surv, sv, family,
                       b, l0i, l0u, w, v, inds){
  rowSums(obs.emp.I(Omega.perturbed, dmats, surv, sv, family, b, l0i, l0u, w, v, inds, TRUE)$Score)
}

.ScoreWrapper <- function(Omega.flat, skeleton, dmats, surv, sv, family,
                          b, l0i, l0u, w, v, inds){
  Omega.perturbed <- relist(Omega.flat, skeleton = skeleton)
  Omega.perturbed$D <- vech2mat(Omega.perturbed$D, sv$q)
  # if(is.not.SPD(Omega.perturbed$D)) Omega.perturbed$D <- pracma::nearest_spd(Omega.perturbed$D)
  .D <<- Omega.perturbed$D
  ScoreOmega(Omega.perturbed, dmats, surv, sv, family,
             b, l0i, l0u, w, v, inds)
}

Xu.Information <- function(Omega, dmats, surv, sv, family,
                           b, l0i, l0u, w, v, inds,
                           method = "central", heps = 1e-4){
  O <- Omega
  O$D <- vech(O$D)
  params <- unlist(as.relistable(O))
  skel <- attr(params, 'skeleton')
  
  numDiff(params, .ScoreWrapper, skeleton = skel, 
          dmats = dmats, surv = surv, sv = sv, family = family,
          b = b, l0i = l0i, l0u = l0u, w = w, v = v, inds = inds,
          method = method, heps = heps)
}

