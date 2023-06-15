## R CMD check results
0 Errors | 0 Warnings | 1 Note.

Note: 
    installed size is  7.6Mb
    sub-directories of 1Mb or more:
      libs   7.2Mb

i.e. gmvjoint.so has size 7.2Mb.

## Test environments 

* Local Ubuntu (R version 3.6.3).
* MacOS (Release, via mac.r-project.org).
* Windows (Release, oldrelease and devel, via win-builder.r-project.org).

R-hub build OK (some NOTEs about "lastMiKTeXException" detritus).

## --dont-test
Some donttest examples have elapsed times > 5s, on my local: 
                        user  system elapsed
   boot.joint        533.070 128.622 389.438
   dynPred            96.688  55.712  40.430
   ROC               105.741  45.484  51.568
   joint              56.362  13.847  43.665
   summary.joint      26.796   4.073  22.530
   residuals.joint    20.078   4.471  15.246
   fitted.joint       18.144   3.987  13.349
   ranef.joint         9.857   1.791   7.932
   anova.joint         9.257   2.127   6.839
   cond.ranefs         7.093   1.733   5.243
   plot.cond.b.joint   6.614   1.310   5.057

