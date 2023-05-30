# gmvjoint version 0.3.0
* Rework to instead use glmmTMB's fits (used for initial conditions) to create required
data matrices instead of a separate step, which frequently created a bottleneck for large
number of responses.
* Added negative binomial distribution.
* Allowed for (possibly time varying) dispersion models (for Gamma/negative binomial/
generalised poisson).

# gmvjoint version 0.2.0
* Dynamic predictions added
* ROC/AUC function added
* Function to generate bootstrapped confidence intervals for AUC added.
* Various S3 methods for joint objects added.
* Relaxation of (quite stringent) model calls.
* Extraction of conditional distribution of random effects.

# gmvjoint version 0.1.0
* First release
