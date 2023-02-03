prepData <- function(data, idVar, timeVar, Yvars, survival, Xvars = NULL){
  if(missing(timeVar)) stop("Please supply timeVar, the longitudinal time variable.")
  if(missing(idVar)) stop("Please supply idVar, the subject identification variable.")
  if(missing(Yvars)) stop("Please supply at least one longitudinal variable")
  if(missing(survival)) stop("Please supply (survival time, event status) in the survival argument.")
  if(!missing(survival) && !class(survival)!='character' && )
  if(missing(data)) stop("Please supply a dataset.")
  
  if (any(c("tbl_df", "tbl") %in% class(data))) data <- as.data.frame(data)
  
  # If Xvars not specified but Yvars then Xvars is everything not specified by other arugments.
  if(is.null(Xvars)) Xvars <- setdiff(names(data), c(idVar, timeVar, Yvars, survival))
  
  data <- data[,c(idVar, timeVar, Xvar, Yvar)]
  
  
}