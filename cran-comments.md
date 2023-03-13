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
NOTE on windows/CRAN submission is on "checking CRAN incoming feasibility", which I can do nothing about.
Release version of windows seems slower than development counterpart for one example (vcov.joint), but both are 
< 10s on my most recent testing.
