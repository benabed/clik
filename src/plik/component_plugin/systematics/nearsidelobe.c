#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void fill_offset_freq_TP(int idreq,double *dreq, int nfreq, double* freqlist,int *mv,int off, error **err);

void nslb_compute(parametric* egl, double *Rq, error **err);

parametric *nslb_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int *mv,m1,m2,f1;
  double dreq[4];
  pfchar name;
  char tp[2];
  
  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 3001*8 + 2*8) + sizeof(int)*egl->nfreq*2,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* (3001*12*12 + 3001*8));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 +2*8);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  if (has_TEB[1]!=0 || has_TEB[2]!=0) {
    fill_offset_freq_TP(4,dreq, egl->nfreq_P,egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
    forwardError(*err,__LINE__,NULL);
  }

  tp[0] = 'T';
  tp[1] = 'P';

  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<2;f1++) {
      sprintf(name,"nslb_epsilon_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
      sprintf(name,"nslb_fwhm_%d_%c",(int)dreq[m1],tp[f1]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);      
    }
  }


  egl->eg_compute = &nslb_compute;
  egl->eg_free = &parametric_simple_payload_free;
    
 return egl;
}


void nslb_compute(parametric* egl, double *Rq, error **err) {
  double *template,*bl,*epsilon,*sigma;
  int *mv;
  double b,v;
  int m1,m2,f1,f2,ell;
  int dreq[4];
  char tp[2];
  pfchar name;
  int m1p,m2p;
  double bp;        

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  tp[0] = 'T';
  tp[1] = 'P';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 3001*8 + 2*8);
  epsilon = egl->payload + sizeof(double)* (3001*12*12 + 3001*8);
  sigma = egl->payload + sizeof(double)* (3001*12*12+ + 3001*8 +8);

  template = egl->payload;
  bl = egl->payload + sizeof(double)* (3001*12*12);

  for(m1=0;m1<4;m1++) {
    for(f1=0;f1<3;f1++) {
      sprintf(name,"nslb_epsilon_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      epsilon[m1*2+f1] = v;
      sprintf(name,"nslb_fwhm_%d_%c",(int)dreq[m1],tp[f1]);
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      sigma[m1*2+f1] = v/sqrt(8*log(2));
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        //if(ell==egl->lmin) {
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        //}
        b = epsilon[mv[m1]]*exp(-ell*(ell+1)*sigma[mv[m1]])/bl[ell*8+mv[m1]] + epsilon[mv[m2]]*exp(-ell*(ell+1)*sigma[mv[m2]])/bl[ell*8+mv[m2]];
        if (m1<egl->nfreq_T*egl->has_TEB[0] && m2>egl->nfreq_T*egl->has_TEB[0]) {
          // TP case, need to symetrise;
          m1p = m1+egl->nfreq_T*egl->has_TEB[0];
          m2p = m2-egl->nfreq_T*egl->has_TEB[0];
          bp = epsilon[mv[m1p]]*exp(-ell*(ell+1)*sigma[mv[m1p]])/bl[ell*8+mv[m1p]] + epsilon[mv[m2p]]*exp(-ell*(ell+1)*sigma[mv[m2p]])/bl[ell*8+mv[m2p]];
          b = .5*(b+bp);
        }
        Rq[IDX_R(egl,ell,m1,m2)] = b*template[ell*12*12+mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
  
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(nslb,nslb_init);
