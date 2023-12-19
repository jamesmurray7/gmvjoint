## R CMD check results

0 Errors | 0 Warnings | 1 Note. 

Note: installed size is 6.9Mb sub-directories of 1Mb or more: libs 6.5Mb

i.e. gmvjoint.so has size 6.5Mb (this NOTE is from both mac and win build test environments). 

1 WARNINGS for mac build test: 

Warning: package `glmmTMB` was built under R version 4.3.1. (from mac build test), Unsure if this is fixable on my end.

## Test environments

-   Local Ubuntu (R version 4.3.2).
-   MacOS (Release, via mac.r-project.org).
-   Windows (Release, oldrelease and devel, via
    win-builder.r-project.org).

(Status: OK for all, besides aforementioned mac WARNING).

## --dont-test

Some `donttest`{...} examples have elapsed times \>5s, on my machine: 

|`fn`|user|system|elapsed|
|:---|:---|:---|:---|
|`boot.joint`|495.930|113.385|364.240|
|`ROC`|116.262|48.197|55.924|
|`dynPred`|103.424|54.159|45.162| 
|`summary.joint`|21.248|4.081|16.584|
|`residuals.joint`|20.035|4.816|14.850|
|`fitted.joint`|18.466|4.151|13.672|
|`ranef.joint`|10.158|1.672|8.097|
|`anova.joint`|8.950|2.096|6.556|
|`plot.cond.b.joint`|8.929|1.763|7.194|
|`cond.ranefs`|8.187|1.770|6.228|
