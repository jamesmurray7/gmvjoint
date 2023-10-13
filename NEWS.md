# gmvjoint version 0.4.0
* Fairly minor updates
* `surv.formula` now properly accepts factor variables. This shortcoming was brought to my attention by J. Piekos (Oxford); so thanks!
* `residuals.joint` functions have been amended to actually work for Cox-Snell (there was an issue with e.g. spline time specifications).
* `plot.residuals.joint` now produces a series of plots with telligible legends.

# gmvjoint version 0.3.0
* Fairly extensive re-writes to workhorse functions `joint` (and underlying `EMUpdate`); some computational efficiency has been gained here 'within' the EM algorithm.
* Removed `CreateDataMatrices` function, in favour of simply taking ones from fits used in initial conditions (via `glmmTMB`); this greatly reduces computation times, particularly for larger models.
* Added a bootstrapping function `boot.joint`: The approximate SEs created by the model are known to be underestimated, and so this bootstrapping function seeks to provide an alternative route; or means of comparison.
* Spruced-up plotting of conditional distributions of the random effects.
* Updates to S3 methods.

# gmvjoint version 0.2.0
* Dynamic predictions added
* ROC/AUC function added
* Function to generate bootstrapped confidence intervals for AUC added.
* Various S3 methods for joint objects added.
* Relaxation of (quite stringent) model calls.
* Extraction of conditional distribution of random effects.

# gmvjoint version 0.1.0
* First release
