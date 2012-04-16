#ifndef _CPADN_
#define _CPADN_
#include "clik_helper.h"

void base_parametric_hdf5_init(hid_t comp_id,char* cur_lkl,int ndet, int** detlist,int *ndef, char ***defkeys, char*** defvalues, int *nvar, char ***varkeys, error **err);
SmicaComp * finalize_parametric_hdf5_init(parametric* p_model,hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err);
#endif
