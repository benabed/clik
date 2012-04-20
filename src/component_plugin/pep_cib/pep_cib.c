#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PEP CIB model
////////////////////////////////////////////////////////////////////////////////////////



void ir_clustered_pep_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double ir_clustered_pep_norm;
  double *A;
  int hfi_freqlist[6] = {100,143,217,353,545,857};
  int nfreqs_hfi = 6;
  int lmax_in = 3000;
  int *ind_freq;
  int ind1, ind2;

  egl = exg;
  A = egl->payload; // Stores input C_l(nu,nu') for all HFI frequencies and ell in [0,3000]
  nfreq = egl->nfreq;
  ind_freq = malloc_err(sizeof(int)*nfreq,err);
  forwardError(*err,__LINE__,);

  ir_clustered_pep_norm = parametric_get_value(egl,"ir_clustered_pep_norm",err);
  forwardError(*err,__LINE__,);


  for (m1=0;m1<nfreq;m1++) {
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (fabs(egl->freqlist[m1]-hfi_freqlist[m2])<1e-6) {
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
        Rq[mell + m1*nfreq+m2] = ir_clustered_pep_norm * A[mell_in + ind1*nfreqs_hfi + ind2];
        Rq[mell + m2*nfreq+m1] = Rq[mell + m1*nfreq+m2];
      }
    }
  }

  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      if (strcmp(egl->varkey[iv],"ir_clustered_pep_norm")==0) {
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = Rq[mell+m1*nfreq+m2]/ir_clustered_pep_norm;
              dRq[mv+mell + m2*nfreq+m1] = Rq[mell+m2*nfreq+m1]/ir_clustered_pep_norm;
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

void ir_clustered_pep_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

parametric *ir_clustered_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int hfi_freqlist[6] = {100,143,217,353,545,857};
  int nfreqs_hfi = 6;
  int m1,m2;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    int ok;
    ok = 0;
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (fabs(egl->freqlist[m1]-hfi_freqlist[m2])<1e-6) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,NULL,egl->freqlist[m1]);
  }

  egl->eg_compute = &ir_clustered_pep_compute;
  egl->eg_free = &ir_clustered_pep_free;

  egl->payload = malloc_err(sizeof(double)*(nfreqs_hfi*nfreqs_hfi*(lmax_in+1)),err); //6 frequencies (HFI), lmax_in=3000
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,rq_clustered_in,nfreqs_hfi*nfreqs_hfi*(lmax_in+1)*sizeof(double));

  parametric_set_default(egl,"ir_clustered_pep_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ir_poisson_pep_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_poisson_pep_norm;
  double *A;
  int hfi_freqlist[6] = {100,143,217,353,545,857};
  int nfreqs_hfi = 6;
  int *ind_freq;
  int ind1, ind2;

  egl = exg;
  A = egl->payload; // Stores input C_l(nu,nu') for all HFI frequencies
  nfreq = egl->nfreq;
  ind_freq = malloc_err(sizeof(int)*nfreq,err);
  forwardError(*err,__LINE__,);

  ir_poisson_pep_norm = parametric_get_value(egl,"ir_poisson_pep_norm",err);
  forwardError(*err,__LINE__,);


  for (m1=0;m1<nfreq;m1++) {
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (fabs(egl->freqlist[m1]-hfi_freqlist[m2])<1e-6) {
        ind_freq[m1]=m2;
      }
    }
  }


  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    lell = ell - egl->lmin;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        ind1 = ind_freq[m1];
        ind2 = ind_freq[m2];
        Rq[mell + m1*nfreq+m2] = ir_poisson_pep_norm * A[ind1*nfreqs_hfi + ind2];
        Rq[mell + m2*nfreq+m1] = Rq[mell + m1*nfreq+m2];
      }
    }
  }

  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      if (strcmp(egl->varkey[iv],"ir_poisson_pep_norm")==0) {
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = Rq[mell+m1*nfreq+m2]/ir_poisson_pep_norm;
              dRq[mv+mell + m2*nfreq+m1] = Rq[mell+m2*nfreq+m1]/ir_poisson_pep_norm;
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

void ir_poisson_pep_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

parametric *ir_poisson_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_poisson_in, error **err) {
  parametric *egl;
  int hfi_freqlist[6] = {100,143,217,353,545,857};
  int nfreqs_hfi = 6;
  int m1,m2;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    int ok;
    ok = 0;
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (fabs(egl->freqlist[m1]-hfi_freqlist[m2])<1e-6) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,NULL,egl->freqlist[m1]);
  }

  egl->eg_compute = &ir_poisson_pep_compute;
  egl->eg_free = &ir_poisson_pep_free;

  egl->payload = malloc_err(sizeof(double)*nfreqs_hfi*nfreqs_hfi,err); //6 frequencies (HFI)
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,rq_poisson_in,nfreqs_hfi*nfreqs_hfi*sizeof(double));
  
  parametric_set_default(egl,"ir_poisson_pep_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_clustered_pep,ir_clustered_pep_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ir_poisson_pep,ir_poisson_pep_init);
