#ifndef _CPADN_
#define _CPADN_
#include "clik_helper.h"

void base_parametric_hdf5_init(hid_t comp_id,char* cur_lkl,int ndet, int** detlist,int *ndef, char ***defkeys, char*** defvalues, int *nvar, char ***varkeys, error **err);
SmicaComp * finalize_parametric_hdf5_init(parametric* p_model,hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err);

#define CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(NAME,INIT_FUNC) \
SmicaComp * clik_smica_comp_##NAME##_init(hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;  \
  SmicaComp *SC;  \
  int lmin,lmax;  \
  int *detlist;  \
  int ndef,nvar;  \
  char **defkeys,**defvalues,**varkeys;  \
  herr_t hstat;  \
  double *template;  \
  int dz;  \
  \
  lmin = ell[0];  \
  lmax = ell[nell-1];  \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);  \
  \
  base_parametric_hdf5_init(comp_id,cur_lkl,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  dz = -1;  \
  template = hdf5_double_datarray(comp_id,cur_lkl,"template",&dz, err);  \
  forwardError(*err,__LINE__,);  \
  p_model = INIT_FUNC(m, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax,template, err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  free(detlist);  \
  if (defkeys[0]!=NULL) {  \
    free(defkeys[0]);  \
    free(defvalues[0]);  \
  }  \
  free(defkeys);  \
  free(defvalues);  \
  \
  if (varkeys[0]!=NULL) {  \
    free(varkeys[0]);  \
  }  \
  free(varkeys);  \
  \
  SC = finalize_parametric_hdf5_init(p_model,comp_id,cur_lkl,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  return SC;  \
}

#define CREATE_PARAMETRIC_FILE_INIT(NAME,INIT_FUNC) \
SmicaComp * clik_smica_comp_##NAME##_init(hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;  \
  SmicaComp *SC;  \
  int lmin,lmax;  \
  int *detlist;  \
  int ndef,nvar;  \
  char **defkeys,**defvalues,**varkeys;  \
  herr_t hstat;  \
  double *template;  \
  int dz;  \
  \
  lmin = ell[0];  \
  lmax = ell[nell-1];  \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);  \
  \
  base_parametric_hdf5_init(comp_id,cur_lkl,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  p_model = INIT_FUNC(m, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax, err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  free(detlist);  \
  if (defkeys[0]!=NULL) {  \
    free(defkeys[0]);  \
    free(defvalues[0]);  \
  }  \
  free(defkeys);  \
  free(defvalues);  \
  \
  if (varkeys[0]!=NULL) {  \
    free(varkeys[0]);  \
  }  \
  free(varkeys);  \
  \
  SC = finalize_parametric_hdf5_init(p_model,comp_id,cur_lkl,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);  \
  forwardError(*err,__LINE__,NULL);  \
    \
  return SC;  \
}

#endif
