#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PEP CIB model
////////////////////////////////////////////////////////////////////////////////////////
void psm_cib_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double psm_cib_norm;
  double *A;
  int hfi_freqlist[4] = {100,143,217,353};
  int nfreqs_hfi = 4;
  int lmax_in = 3000;
  int *ind_freq;
  int ind1, ind2;

  egl = exg;
  A = egl->payload; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = malloc_err(sizeof(int)*nfreq,err);
  forwardError(*err,__LINE__,);

  psm_cib_norm = parametric_get_value(egl,"psm_cib_norm",err);
  forwardError(*err,__LINE__,);


  for (m1=0;m1<nfreq;m1++) {
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (egl->freqlist[m1] == hfi_freqlist[m2]) {
        ind_freq[m1]=m2;
      }
    }
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    mell_in=ell*nfreqs_hfi*nfreqs_hfi;
    lell = ell - egl->lmin;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[mell + m1*nfreq+m2] = psm_cib_norm *1e12* A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[mell + m2*nfreq+m1] = Rq[mell + m1*nfreq+m2];
      }
    }
  }

  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      if (strcmp(egl->varkey[iv],"psm_cib_norm")==0) {
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = Rq[mell+m1*nfreq+m2]/psm_cib_norm;
              dRq[mv+mell + m2*nfreq+m1] = Rq[mell+m2*nfreq+m1]/psm_cib_norm;
            }
          }
        }
        continue;
      }
      // error return
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }

  free(ind_freq);
  return;
}

void psm_cib_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

parametric *psm_cib_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int hfi_freqlist[4] = {100,143,217,353};
  int nfreqs_hfi = 4;
  int m1,m2;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    int ok;
    ok = 0;
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (egl->freqlist[m1]==hfi_freqlist[m2]) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,NULL,egl->freqlist[m1]);
  }

  egl->eg_compute = &psm_cib_compute;
  egl->eg_free = &psm_cib_free;

  egl->payload = malloc_err(sizeof(double)*(nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),err); //6 frequencies (HFI), lmax_in=3000
  forwardError(*err,__LINE__,NULL);
  
  memcpy(egl->payload,rq_clustered_in,nfreqs_hfi*nfreqs_hfi*(lmax_in+1)*sizeof(double));

  parametric_set_default(egl,"psm_cib_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}


SmicaComp * clik_smica_comp_psm_cib_init(hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  parametric* p_model;
  SmicaComp *SC;
  int lmin,lmax;
  int *detlist;
  int ndef,nvar;
  char **defkeys,**defvalues,**varkeys;
  herr_t hstat;
  double *template;
  int dz;

  lmin = ell[0];
  lmax = ell[nell-1];
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);

  base_parametric_hdf5_init(comp_id,cur_lkl,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);
  forwardError(*err,__LINE__,NULL);
  
  dz = -1;
  template = hdf5_double_datarray(comp_id,cur_lkl,"template",&dz, err);
  forwardError(*err,__LINE__,);

  p_model = psm_cib_init(m, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax,template, err);
  forwardError(*err,__LINE__,NULL);
  
  free(detlist);
  if (defkeys[0]!=NULL) {
    free(defkeys[0]);
    free(defvalues[0]);
  }
  free(defkeys);
  free(defvalues);

  if (varkeys[0]!=NULL) {
    free(varkeys[0]);
  }
  free(varkeys);

  SC = finalize_parametric_hdf5_init(p_model,comp_id,cur_lkl,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);
  forwardError(*err,__LINE__,NULL);
  
  return SC;
}