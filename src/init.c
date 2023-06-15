#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gmvjoint_appxE_Gammasigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_appxE_GenPoissigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_appxE_NegBinsigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_dmvn_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_dmvt_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Egammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_GP1_pmf_scalar(SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Hbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Hgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_joint_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_joint_density_ddb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_lambda_hat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_lambda_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_logfti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_make_eta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_make_tau(SEXP, SEXP);
extern SEXP _gmvjoint_metropolis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_S_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Sbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Sgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_sigma2_Gaussian_update(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_vech2mat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gmvjoint_appxE_Gammasigma",       (DL_FUNC) &_gmvjoint_appxE_Gammasigma,        7},
    {"_gmvjoint_appxE_GenPoissigma",     (DL_FUNC) &_gmvjoint_appxE_GenPoissigma,      7},
    {"_gmvjoint_appxE_NegBinsigma",      (DL_FUNC) &_gmvjoint_appxE_NegBinsigma,       7},
    {"_gmvjoint_dmvn_fast",              (DL_FUNC) &_gmvjoint_dmvn_fast,               4},
    {"_gmvjoint_dmvt_fast",              (DL_FUNC) &_gmvjoint_dmvt_fast,               5},
    {"_gmvjoint_Egammazeta",             (DL_FUNC) &_gmvjoint_Egammazeta,             13},
    {"_gmvjoint_GP1_pmf_scalar",         (DL_FUNC) &_gmvjoint_GP1_pmf_scalar,          3},
    {"_gmvjoint_Hbeta",                  (DL_FUNC) &_gmvjoint_Hbeta,                  14},
    {"_gmvjoint_Hgammazeta",             (DL_FUNC) &_gmvjoint_Hgammazeta,             14},
    {"_gmvjoint_joint_density",          (DL_FUNC) &_gmvjoint_joint_density,          21},
    {"_gmvjoint_joint_density_ddb",      (DL_FUNC) &_gmvjoint_joint_density_ddb,      21},
    {"_gmvjoint_lambda_hat",             (DL_FUNC) &_gmvjoint_lambda_hat,              9},
    {"_gmvjoint_lambda_update",          (DL_FUNC) &_gmvjoint_lambda_update,          10},
    {"_gmvjoint_logfti",                 (DL_FUNC) &_gmvjoint_logfti,                 10},
    {"_gmvjoint_make_eta",               (DL_FUNC) &_gmvjoint_make_eta,                6},
    {"_gmvjoint_make_tau",               (DL_FUNC) &_gmvjoint_make_tau,                2},
    {"_gmvjoint_metropolis",             (DL_FUNC) &_gmvjoint_metropolis,             23},
    {"_gmvjoint_S_",                     (DL_FUNC) &_gmvjoint_S_,                      4},
    {"_gmvjoint_Sbeta",                  (DL_FUNC) &_gmvjoint_Sbeta,                  14},
    {"_gmvjoint_Sgammazeta",             (DL_FUNC) &_gmvjoint_Sgammazeta,             14},
    {"_gmvjoint_sigma2_Gaussian_update", (DL_FUNC) &_gmvjoint_sigma2_Gaussian_update,  5},
    {"_gmvjoint_vech2mat",               (DL_FUNC) &_gmvjoint_vech2mat,                2},
    {NULL, NULL, 0}
};

void R_init_gmvjoint(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
