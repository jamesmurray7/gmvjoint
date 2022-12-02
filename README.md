# `GMVJM`

## What is `GMVJM`?

* "**G**eneralised"
* **M**ulti**v**ariate
* **J**oint **M**odels.


GMVJM allows the user to fit joint models of survival and multivariate longitudinal data, where the 
longitudinal sub-models are specified by generalised linear mixed models (GLMMs). The joint models 
are fit via maximum likelihood using an approximate EM algorithm first proposed by Bernhardt *et
al*. (2015). The GLMMs are specified using the same syntax as for package `glmmTMB` (Brooks *et
al*., 2017). The joint models themselves are then the  flexible extensions to those in e.g.
Wulfoshn and Tsiatis (1997). The user is able to simulate data under many different response
types.

## To-do list
The package in current incantation is relatively skeletal, as such not a lot of post-hoc
anaylses on fitted joint models is possible. As such, an immediate to-do list currently looks like

* Dynamic predictions: Note that code to do this already exists at https://github.com/jamesmurray7/GLM/blob/main/Multi-test/DynamicPredictions.R, 
which just needs to be ported over.
* Usual modelling generics: Fitted, residuals (probably just longitudinal) and plots. 

No promises are made w.r.t timescale of these being implemented: Currently I am a PhD
student and little cache is awarded for production or maintenance of R packages!

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

