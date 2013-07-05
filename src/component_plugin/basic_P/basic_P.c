#include "clik_parametric.h"
#include "clik_parametric_addon.h"

typedef struct {
  int n;
  int *m1,*m2;
  double *nrm;
  char **nrm_names;
} pw_XX_payload;

void pw_TE_compute(parametric* exg, double *Rq, error **err);

void pw_XX_free(void **PP) {
  pw_XX_payload *P;

  P = *PP;
  free(P->nrm);
  free(P->m1);
  free(P->m2);
  free(P->nrm_names[0]);
  free(P->nrm_names);
  free(P);
  *PP = NULL;
}

parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int hk;
  char **chn, *bf,*cbf,*pch;
  char *channeln;
  int tl, mm, lbf,l,delta;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &pw_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = malloc_err(sizeof(pw_XX_payload),err);
  forwardError(*err,__LINE__,NULL);
    
  payload = egl->payload;
  payload->n = egl->nfreq_T * egl->nfreq_P;
  mx = egl->nfreq_T + egl->nfreq_P;

  payload->nrm = malloc_err(sizeof(double)*mx,err);
  forwardError(*err,__LINE__,NULL);
  payload->m1 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  payload->m2 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  
  mx = 0;
  for(i=0;i<egl->nfreq_T;i++) {
    for(j=0;j<egl->nfreq_P;j++) {
      payload->m1[mx] = i;
      payload->m2[mx] = j + egl->nfreq_T;
      mx++;
    }
  }

  mx = egl->nfreq_T + egl->nfreq_P;

  parametric_set_default(egl,"pw_TE_index",-1,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pw_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);
  
  bf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  cbf = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);
  chn = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<egl->nfreq_T;i++) {
    sprintf(cbf,"%s,%d_T",cbf,i);
  }
  for(i=0;i<egl->nfreq_P;i++) {
    sprintf(cbf,"%s,%d_E",cbf,i);
  }

  hk = cdic_key_index(egl->pf,"pw_TE_channel_names",err);
  forwardError(*err,__LINE__,NULL);
  
  channeln = cbf;

  
  if (hk!=-1) {
    channeln = cdic_get(egl->pf,"pw_TE_channel_names",NULL,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  tl = strlen(channeln);
  mm = 0;
  lbf = 0;
  l=0;
  while(l<tl) {
    chn[mm] = bf + lbf;
    while(channeln[l]==' ' || channeln[l]==',') {
      l++;
    }
    while(l<tl && channeln[l]!=',') {
      bf[lbf] = channeln[l];
      l++;
      lbf++;
    }
    bf[lbf] = '\0';
    lbf++;
    //_DEBUGHERE_("%d |%s|",mm,chn[mm]);
    mm++;
  }
  testErrorRet(mm!=mx,-97664,"bad pw_TE_channel_names value",*err,__LINE__,NULL);

  pch = malloc_err(sizeof(char)*mx*256,err);
  forwardError(*err,__LINE__,NULL);  

  payload->nrm_names = malloc_err(sizeof(char*)*mx,err);
  forwardError(*err,__LINE__,NULL);  

  delta = 0;

  for(i=0;i<mx;i++) {
    payload->nrm_names[i] = pch + delta;
    sprintf(payload->nrm_names[i],"pw_TE_A_%s",chn[i]);
    parametric_set_default(egl,payload->nrm_names[i],1.,err);
    forwardError(*err,__LINE__,NULL); 
   
    delta = strlen(payload->nrm_names[i]) +1;
    pch = payload->nrm_names[i];
  }      

  free(bf);
  free(cbf);
  return egl;
}

void pw_TE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v;

  payload = egl->payload;
  mx = egl->nfreq_T + egl->nfreq_P;
  

  for(i=0;i<mx;i++) {
    payload->nrm[i] = parametric_get_value(egl,payload->nrm_names[i],err);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"pw_TE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"pw_TE_index",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

CREATE_PARAMETRIC_POL_FILE_INIT(pw_TE,pw_TE_init);
