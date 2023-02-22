#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gmvjoint_dmvnrm_arma_fast(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_dmvt_arma_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Egammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_GP1_pmf_scalar(SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Hbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Hgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_joint_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_joint_density_ddb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_joint_density_sdb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_lambdaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_ll_Gamma(SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_ll_genpois(SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_logfti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_maketau(SEXP, SEXP);
extern SEXP _gmvjoint_phi_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_S_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Sbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_Sgammazeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_vare_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gmvjoint_vech2mat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gmvjoint_dmvnrm_arma_fast",  (DL_FUNC) &_gmvjoint_dmvnrm_arma_fast,   4},
    {"_gmvjoint_dmvt_arma_fast",    (DL_FUNC) &_gmvjoint_dmvt_arma_fast,     5},
    {"_gmvjoint_Egammazeta",        (DL_FUNC) &_gmvjoint_Egammazeta,        13},
    {"_gmvjoint_GP1_pmf_scalar",    (DL_FUNC) &_gmvjoint_GP1_pmf_scalar,     3},
    {"_gmvjoint_Hbeta",             (DL_FUNC) &_gmvjoint_Hbeta,             13},
    {"_gmvjoint_Hgammazeta",        (DL_FUNC) &_gmvjoint_Hgammazeta,        15},
    {"_gmvjoint_joint_density",     (DL_FUNC) &_gmvjoint_joint_density,     20},
    {"_gmvjoint_joint_density_ddb", (DL_FUNC) &_gmvjoint_joint_density_ddb, 20},
    {"_gmvjoint_joint_density_sdb", (DL_FUNC) &_gmvjoint_joint_density_sdb, 21},
    {"_gmvjoint_lambdaUpdate",      (DL_FUNC) &_gmvjoint_lambdaUpdate,      10},
    {"_gmvjoint_ll_Gamma",          (DL_FUNC) &_gmvjoint_ll_Gamma,           3},
    {"_gmvjoint_ll_genpois",        (DL_FUNC) &_gmvjoint_ll_genpois,         3},
    {"_gmvjoint_logfti",            (DL_FUNC) &_gmvjoint_logfti,            10},
    {"_gmvjoint_maketau",           (DL_FUNC) &_gmvjoint_maketau,            2},
    {"_gmvjoint_phi_update",        (DL_FUNC) &_gmvjoint_phi_update,         9},
    {"_gmvjoint_S_",                (DL_FUNC) &_gmvjoint_S_,                 4},
    {"_gmvjoint_Sbeta",             (DL_FUNC) &_gmvjoint_Sbeta,             13},
    {"_gmvjoint_Sgammazeta",        (DL_FUNC) &_gmvjoint_Sgammazeta,        14},
    {"_gmvjoint_vare_update",       (DL_FUNC) &_gmvjoint_vare_update,        8},
    {"_gmvjoint_vech2mat",          (DL_FUNC) &_gmvjoint_vech2mat,           2},
    {NULL, NULL, 0}
};

void R_init_gmvjoint(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
