## R CMD check results
0 Errors | 0 Warnings | 1 Note.

Note: 
   installed size is  8.7Mb
   sub-directories of 1Mb or more:
     libs   8.3Mb

i.e. gmvjoint.so has size 8.3Mb.

## Test environments 

* Local Ubuntu (R version 3.6.3).
* MacOS (Release, via mac.r-project.org).
* Windows (Release and devel, via win-builde.r-project.org).

R-hub build OK.
Note upon installation one may be warned about inconsistencies across TMB/Matrix versions, but this
is not within my control.
Commented warnings about CXX_STD=CXX11 in Makevars(/.win) which resulted in a NOTE. Hopefully this doesn't result
in deterioration of performance since I thought RcppArmadillo required this!