#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* SO I DONT FORGET: 
   R code to generate -->
   tools::package_native_routine_registration_skeleton(dir = '.', con = './src/init.c', character_only = F)
*/

/* .Call calls */
extern SEXP _GMVJM_GP1_pmf_scalar(SEXP, SEXP, SEXP);
extern SEXP _GMVJM_Hbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_Hess_eta_genpois(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_Hgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_joint_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_joint_density_ddb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_joint_density_sdb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_lambdaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_ll_genpois(SEXP, SEXP, SEXP);
extern SEXP _GMVJM_logfti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_phi_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_Sbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_Sgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GMVJM_vare_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GMVJM_GP1_pmf_scalar",    (DL_FUNC) &_GMVJM_GP1_pmf_scalar,     3},
    {"_GMVJM_Hbeta",             (DL_FUNC) &_GMVJM_Hbeta,              9},
    {"_GMVJM_Hess_eta_genpois",  (DL_FUNC) &_GMVJM_Hess_eta_genpois,   4},
    {"_GMVJM_Hgammazeta",        (DL_FUNC) &_GMVJM_Hgammazeta,        15},
    {"_GMVJM_joint_density",     (DL_FUNC) &_GMVJM_joint_density,     20},
    {"_GMVJM_joint_density_ddb", (DL_FUNC) &_GMVJM_joint_density_ddb, 20},
    {"_GMVJM_joint_density_sdb", (DL_FUNC) &_GMVJM_joint_density_sdb, 21},
    {"_GMVJM_lambdaUpdate",      (DL_FUNC) &_GMVJM_lambdaUpdate,      13},
    {"_GMVJM_ll_genpois",        (DL_FUNC) &_GMVJM_ll_genpois,         3},
    {"_GMVJM_logfti",            (DL_FUNC) &_GMVJM_logfti,            10},
    {"_GMVJM_phi_update",        (DL_FUNC) &_GMVJM_phi_update,         9},
    {"_GMVJM_Sbeta",             (DL_FUNC) &_GMVJM_Sbeta,              9},
    {"_GMVJM_Sgammazeta",        (DL_FUNC) &_GMVJM_Sgammazeta,        15},
    {"_GMVJM_vare_update",       (DL_FUNC) &_GMVJM_vare_update,        8},
    {NULL, NULL, 0}
};

void R_init_GMVJM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
