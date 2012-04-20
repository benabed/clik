#include "clik_parametric.h"
#include "clik_parametric_addon.h"

////////////////////////////////////////////////////////////////////////////////////////
// PSM CIB rigid model
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


////////////////////////////////////////////////////////////////////////////////////////
// PSM CIB free model
////////////////////////////////////////////////////////////////////////////////////////
void psm_cib_free_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv,lell,mell_in;
  double psm_cib_norm;
  double *A,*N;
  int lmax_in = 3000;
  pfchar name;

  egl = exg;
  A = egl->payload; // Stores input C_l
  nfreq = egl->nfreq;
  N = malloc_err(sizeof(double)*nfreq,err);
  forwardError(*err,__LINE__,);

  for(m1=0;m1<nfreq;m1++) {    
    sprintf(name,"A_%d",egl->freqlist[m1]);
    N[m1] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,);
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    lell = ell - egl->lmin;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[mell + m1*nfreq+m2] = N[m1]*N[m2] * A[ell];
        Rq[mell + m2*nfreq+m1] = Rq[mell + m1*nfreq+m2];
      }
    }
  }

  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      int m,done;
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      done = 0;
      for(m=0;m<nfreq;m++) {
        sprintf(name,"A_%d",egl->freqlist[m]);
        if (strcmp(egl->varkey[iv],name)==0) {
          memset(dRq,0,sizeof(double)*(egl->lmax-egl->lmin+1)*nfreq*nfreq);
          for (ell=egl->lmin;ell<=egl->lmax;ell++) {
            mell = (ell-egl->lmin)*nfreq*nfreq;
            for (m2=0;m2<nfreq;m2++) {
              if (m2==m) {
                dRq[mv+mell+m*nfreq+m] = 2*Rq[mell+m*nfreq+m]/N[m];
                continue;    
              }
              dRq[mv+mell + m*nfreq+m2] = Rq[mell+m*nfreq+m2]/N[m];
              dRq[mv+mell + m2*nfreq+m] = dRq[mv+mell + m*nfreq+m2];
            }
          }
          done=1;
          break;
        }
      }
      if (done==1) {
        continue;
      }
      // error return
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }

  free(N);
  return;
}


parametric *psm_cib_free_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;
  int lmax_in=3000;
  int m1;
  pfchar name;

  testErrorRetVA(lmax>lmax_in,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,NULL,lmax,lmax_in);
   
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
   
  egl->eg_compute = &psm_cib_free_compute;
  egl->eg_free = &psm_cib_free;
   
  egl->payload = malloc_err(sizeof(double)*((lmax_in+1)),err); //6 frequencies (HFI), lmax_in=3000
  forwardError(*err,__LINE__,NULL);
   
  memcpy(egl->payload,rq_clustered_in,(lmax_in+1)*sizeof(double));
   
  for(m1=0;m1<egl->nfreq;m1++) {
    sprintf(name,"A_%d",egl->freqlist[m1]);
    pt = name;

    parametric_set_default(egl,pt,1,err);
    forwardError(*err,__LINE__,NULL);
  }
   
  return egl;
}


CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_cib,psm_cib_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(psm_cib_free,psm_cib_free_init);
