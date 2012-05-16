#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM diffuse emission rigid model, by detector/detset
////////////////////////////////////////////////////////////////////////////////////////
void psm_diffuse_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,ndet,iv,mv,lell,mell_in;
  double psm_diffuse_bydet_norm;
  double *A;
  int ndets_hfi = 19;
  int lmax_in = 3000;
  template_payload *payload;

  payload = egl->payload;
  A = payload->template; // Stores input C_l(det,det') for det,det' in all HFI detsets/detectors and ell in [0,3000]
  ndet = egl->ndet;
  ind_freq = payload->ind_freq;

  psm_diffuse_bydet_norm = parametric_get_value(egl,"psm_diffuse_bydet_norm",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell_in=ell*ndets_hfi*ndets_hfi;
    for (m1=0;m1<ndet;m1++) {
      for (m2=m1;m2<ndet;m2++) {
        ind1 = ind_freq[m1]; // Think "ind_det", just recycling payload template structure
        ind2 = ind_freq[m2]; // Think "ind_det", just recycling payload template structure
        Rq[IDX_R(egl,ell,m1,m2)] = psm_diffuse_bydet_norm *1e12* A[mell_in + ind1*ndets_hfi + ind2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

parametric *psm_diffuse_bydet_init(int ndet, char **detnames, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  char hfi_detnames[19][7] = {"100_ds1","100_ds2","143_ds1","143_ds2","143_5","143_6","143_7","217_1","217_2","217_3","217_4","217_ds1","217_ds2","353_1","353_2","353_ds1","353_7","353_8"};
  int ndets_hfi = 19;
  int m1,m2;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  //egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl = parametric_bydet_init(ndet, detnames, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  parametric_template_payload_bydet_init(egl,template, (ndets_hfi*ndets_hfi*(lmax_in+1)),hfi_detnames,ndets_hfi,err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &psm_diffuse_bydet_compute;
  egl->eg_free = &parametric_template_payload_free;

  parametric_set_default(egl,"psm_diffuse_bydet_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"psm_diffuse_bydet_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_diffuse_bydet,psm_diffuse_bydet_init);
