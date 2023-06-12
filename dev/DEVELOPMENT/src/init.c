#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DEVgmvjoint_appxE_Gammasigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_appxE_GenPoissigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_appxE_NegBinsigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_dmvn_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_dmvt_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_Egammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_GP1_pmf_scalar(SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_Hbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_Hgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_joint_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_joint_density_ddb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_lambda_hat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_lambda_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_logfti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_make_eta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_make_tau(SEXP, SEXP);
extern SEXP _DEVgmvjoint_metropolis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_S_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_Sbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_Sgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_sigma2_Gaussian_update(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DEVgmvjoint_vech2mat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_DEVgmvjoint_appxE_Gammasigma",       (DL_FUNC) &_DEVgmvjoint_appxE_Gammasigma,        7},
    {"_DEVgmvjoint_appxE_GenPoissigma",     (DL_FUNC) &_DEVgmvjoint_appxE_GenPoissigma,      7},
    {"_DEVgmvjoint_appxE_NegBinsigma",      (DL_FUNC) &_DEVgmvjoint_appxE_NegBinsigma,       7},
    {"_DEVgmvjoint_dmvn_fast",              (DL_FUNC) &_DEVgmvjoint_dmvn_fast,               4},
    {"_DEVgmvjoint_dmvt_fast",              (DL_FUNC) &_DEVgmvjoint_dmvt_fast,               5},
    {"_DEVgmvjoint_Egammazeta",             (DL_FUNC) &_DEVgmvjoint_Egammazeta,             13},
    {"_DEVgmvjoint_GP1_pmf_scalar",         (DL_FUNC) &_DEVgmvjoint_GP1_pmf_scalar,          3},
    {"_DEVgmvjoint_Hbeta",                  (DL_FUNC) &_DEVgmvjoint_Hbeta,                  14},
    {"_DEVgmvjoint_Hgammazeta",             (DL_FUNC) &_DEVgmvjoint_Hgammazeta,             14},
    {"_DEVgmvjoint_joint_density",          (DL_FUNC) &_DEVgmvjoint_joint_density,          21},
    {"_DEVgmvjoint_joint_density_ddb",      (DL_FUNC) &_DEVgmvjoint_joint_density_ddb,      21},
    {"_DEVgmvjoint_lambda_hat",             (DL_FUNC) &_DEVgmvjoint_lambda_hat,              9},
    {"_DEVgmvjoint_lambda_update",          (DL_FUNC) &_DEVgmvjoint_lambda_update,          10},
    {"_DEVgmvjoint_logfti",                 (DL_FUNC) &_DEVgmvjoint_logfti,                 10},
    {"_DEVgmvjoint_make_eta",               (DL_FUNC) &_DEVgmvjoint_make_eta,                6},
    {"_DEVgmvjoint_make_tau",               (DL_FUNC) &_DEVgmvjoint_make_tau,                2},
    {"_DEVgmvjoint_metropolis",             (DL_FUNC) &_DEVgmvjoint_metropolis,             23},
    {"_DEVgmvjoint_S_",                     (DL_FUNC) &_DEVgmvjoint_S_,                      4},
    {"_DEVgmvjoint_Sbeta",                  (DL_FUNC) &_DEVgmvjoint_Sbeta,                  14},
    {"_DEVgmvjoint_Sgammazeta",             (DL_FUNC) &_DEVgmvjoint_Sgammazeta,             14},
    {"_DEVgmvjoint_sigma2_Gaussian_update", (DL_FUNC) &_DEVgmvjoint_sigma2_Gaussian_update,  5},
    {"_DEVgmvjoint_vech2mat",               (DL_FUNC) &_DEVgmvjoint_vech2mat,                2},
    {NULL, NULL, 0}
};

void R_init_DEVgmvjoint(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
