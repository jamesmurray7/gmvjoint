## R CMD check results
0 Errors | 0 Warnings | 1 Note.

Note: 
   installed size is  7.7Mb
   sub-directories of 1Mb or more:
     libs   7.3Mb

i.e. gmvjoint.so has size 7.3Mb.

## Test environments 

* Local Ubuntu (R version 3.6.3).
* MacOS (Release, via mac.r-project.org).
* Windows (Release, oldrelease and devel, via win-builder.r-project.org).

R-hub build OK (some NOTEs about "lastMiKTeXException" detritus).

## Other Notes
Removed usage to Matrix package, replacing with bespoke code or using another,
which should alleviate deprecation warnings. They may persist for some time due to the 
Matrix package being used in glmmTMB.
Commented warnings about CXX_STD=CXX11 in Makevars(/.win) which resulted in a NOTE. 
Tried some fixes to lambdaUpdate in src/funs.cpp, which was causing a clang issue. Hopefully it remedies.

## --dont-test
Some donttest examples have elapsed times > 5s, on my local: 
                      user system elapsed
   bootAUCdiff     181.121  2.419 178.822
   ROC             123.861 52.477  60.891
   dynPred         111.482 56.890  52.063
   bootAUC         159.765  1.109 158.826
   joint            30.694 12.269  20.309
   residuals.joint  19.426  4.590  15.021
   fitted.joint     19.165  4.317  14.465
   summary.joint    13.257  2.827  10.251
   anova.joint      11.309  2.276   8.798
   cond.ranefs       9.872  2.389   7.308
   ranef.joint      10.523  1.728   8.595

