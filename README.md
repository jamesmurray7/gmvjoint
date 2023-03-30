# `gmvjoint`

## What is `gmvjoint`?

* "**g**eneralised"
* **m**ulti**v**ariate
* **Joint** models.


`gmvjoint` allows the user to fit joint models of survival and multivariate longitudinal data, where the 
longitudinal sub-models are specified by generalised linear mixed models (GLMMs). The joint models 
are fit via maximum likelihood using an approximate EM algorithm first proposed by Bernhardt *et
al*. (2015). The GLMMs are specified using the same syntax as for package `glmmTMB` (Brooks *et
al*., 2017). The joint models themselves are then the  flexible extensions to those in e.g.
Wulfsohn and Tsiatis (1997). The user is able to simulate data under many different response
types.

Currently, five families can be fit: Gaussian; Poisson; binomial; Gamma and generalised Poisson. 

## Installation
You can install the latest 'official' release from CRAN in the usual way: 
```r
install.packages('gmvjoint')
```

or the latest development version using `devtools`: 
```r
devtools::install_github('jamesmurray7/gmvjoint')
``` 

## Example
To fit a joint model, we first need to specify the longitudinal and survival sub-models. 

The longitudinal sub-model **must** be a list which contains the specification of the longitudinal process along with its random effects structure 
in the same syntax as a [glmmTMB](https://cran.r-project.org/package=glmmTMB) model (which itself is the same as the widely-used `lme4`). 
As an example, suppose we want to fit a trivariate model on the oft-used PBC data, with a linear time-drug interaction term on albumin, a spline term on
(logged) serum bilirubin and a linear fit on spiders, we specify
```r
data(PBC)
PBC$serBilir <- log(PBC$serBilir)
PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
                              'serBilir', 'albumin', 'spiders'))
PBC <- na.omit(PBC) 
long.formulas <- list(
  albumin ~ drug * time + (1 + time|id),
  serBilir ~ drug * splines::ns(time, 3) + (1 + splines::ns(time, 3)|id),
  spiders ~ drug * time + (1|id)
)
```
where we note interactions and spline-time fits are possible. Currently, transformations on variables (e.g. `log(Y)`) must be done *before* this setup of formulae. 

The survival sub-model must be set-up using `Surv()` from the [survival](https://cran.r-project.org/package=survival) package e.g.
```r
surv.formula <- Surv(survtime, status) ~ drug
```
Currently interaction terms in the survival sub-model specification are unsupported. 

Now we can do the joint model call through the main workhorse function `joint`. This notably take a *list* of family arguments which **must** match-up in the desired order as the longitudinal process
list. We call our `fit` via
```r
fit <- joint(long.formulas = long.formulas, surv.formula = surv.formula, data = PBC, 
             family = list("gaussian", "gaussian", "binomial"))
summary(fit)
```
where extra control arguments are documented in `?joint`. Numerous S3 methods exist for the class of object `joint` creates: `summary()`, `logLik()`, `fixef()`, `ranef()`, `fitted()` and `resid()`. 

We bridge from a set of joint model parameter estimates to a prognostic one by dynamic predictions `dynPred`. We can assess discriminatory capabilities of the `joint()` model fit by the `ROC` function, too.

## To-do list
Currently the largest limitation exists with the relatively strict data structure necessary and the corresponding calls to the `joint` function. The below lists these (known) limitations and plans for relaxing.

* Longitudinal information: The longitudinal time argument must be named `time` and the subject identifier (which we 'split' random effects by) `id`. Unsure if I will ever change these; I think a little more user pre-processing is no bad thing, when alternative would be a more crowded call to `joint`, which I wouldn't be a fan of.
* Misc.: data must be balanced (i.e. no `NA` values); this will be fixed in a future update. For now I don't think this is the biggest issue, and recommend using `na.omit` for example. Additionally, the id variable __must__ increment by no more than one. That is, `data$id=1,1,1,2,2,2,3,3,3` is fine, but `data$id=1,1,1,1,3,3,3,4,4` is not. This is due to how data matrices are created internally and will be fixed in the future. 
* Something of a bottleneck is creation of data matrices for large K. Some promising results in e.g. dev/CreateDataMatrices_gTMB.R. May be incorporated in future!
* Pretty-ing of progress bars using cli package.

## References

Bernhardt PW, Zhang D and Wang HJ. A fast EM Algorithm for Fitting Joint Models of a Binary 
Response to Multiple Longitudinal Covariates Subject to Detection Limits. 
*Computational Statistics and Data Analysis* 2015; **85**; 37--53

Mollie E. Brooks, Kasper Kristensen, Koen J. van Benthem, Arni Magnusson, Casper W. Berg, Anders
Nielsen, Hans J. Skaug, Martin Maechler and Benjamin M. Bolker (2017). glmmTMB Balances Speed and
Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling. 
*The R Journal*, **9(2)**, 378-400.

Murray, J and Philipson P. A fast approximate EM algorithm for joint models of survival and
multivariate longitudinal data. *Computational Statistics and Data Analysis* 2022

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
measured with error. *Biometrics.* 1997; **53(1)**, 330-339.

