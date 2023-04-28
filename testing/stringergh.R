stringer.gh <- function(data, idcol = 'id', timecol = 'time', 
                        longit.pc = .25){
  data$.temp.id.col <- data[,idcol]; data$.temp.time.col <- data[,timecol]
  n <- length(unique(data$.temp.id.col))
  if(longit.pc > 1 | longit.pc < 0) stop("longit.pc must be in [0,1].")
  r.lp <- quantile(with(data, tapply(.temp.time.col, .temp.id.col, length)),
                  probs = longit.pc)
  return(ceiling(1.5 * log(n, base = r.lp) - 2))
}