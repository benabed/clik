#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def, error **err);

void mul0_compute(parametric *egl, double *rq, error **err);

parametric *mul0_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &mul0_compute;
  
  parametric_set_default(egl,"mul0_value",100,err); 
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

void mul0_compute(parametric *egl, double *Rq, error **err) {
  double l_pivot,delta_l,v;
  int m1,m2,ell;
  double dip_nrm;
  int lmin,lmax;

  v = parametric_get_value(egl,"mul0_value",err);
  forwardError(*err,__LINE__,);
  lmin =  egl->lmin;
  lmax =  egl->lmax;
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=0;m2<egl->nfreq;m2++) {
      for(ell=lmin;ell<=lmax;ell++) {
        Rq[IDX_R(egl,ell,m1,m2)] = v; 
      }
    }
  }
}
CREATE_PARAMETRIC_FILE_INIT(mul0,mul0_init);

void beamnl_compute(parametric *egl, double *rq, error **err);

parametric *beamnl_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  double *Ip;
  int *mv;
  int nfreq_template,lmax_template;
  pfchar name;
  int nmod;
  double *dreq;
  int i,m1,p;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &beamnl_compute;
  egl->eg_free = &parametric_simple_payload_free;

  nfreq_template = parametric_get_value(egl,"beamnl_nfreq_template",err);
  forwardError(*err,__LINE__,NULL);

  lmax_template = parametric_get_value(egl,"beamnl_lmax_template",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"beamnl_nmode",3 ,err); 
  forwardError(*err,__LINE__,NULL);
  nmod = parametric_get_value(egl,"beamnl_nmode",err);
  forwardError(*err,__LINE__,NULL);
  nmod+=1;

  parametric_set_default(egl,"beamnl_strict",1 ,err); 
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*nmod + sizeof(double)*(nmod*2)+sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(egl->payload,template,sizeof(double)* ((lmax_template+1)*(nfreq_template)));
  
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*nmod + sizeof(double)*(nmod*2);
  
  dreq = malloc_err(sizeof(double)*nfreq_template,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<nfreq_template;i++) {
    sprintf(name,"beamnl_freq_%d",i);
    dreq[i] = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,NULL);
  }  

  fill_offset_freq(nfreq_template,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(p=1;p<nmod;p++) {
      sprintf(name,"beamnl_I_%d_%d",(int)egl->freqlist[m1],p);
      parametric_set_default(egl,name,1,err); 
      forwardError(*err,__LINE__,NULL);
    }
  }

  free(dreq);

  return egl;
}

void beamnl_compute(parametric *egl, double *Rq, error **err) {
  double l_pivot,delta_l,v;
  int m1,m2,ell,p,p1,p2;
  double dip_nrm;
  int lmin,lmax;
  double *template,*Ip,*tp1,*tp2;
  int *mv;
  int nfreq_template,lmax_template;
  pfchar name;
  int nmod;
  int strict;
  double t1,t2,lt1,lt2;


  nfreq_template = parametric_get_value(egl,"beamnl_nfreq_template",err);
  forwardError(*err,__LINE__,);

  lmax_template = parametric_get_value(egl,"beamnl_lmax_template",err);
  forwardError(*err,__LINE__,);
  
  nmod = parametric_get_value(egl,"beamnl_nmode",err);
  forwardError(*err,__LINE__,);
  nmod+=1;

  strict = parametric_get_value(egl,"beamnl_strict",err);
  forwardError(*err,__LINE__,);
  
  template = egl->payload;
  Ip = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template));
  tp1 = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*nmod;
  tp2 = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*nmod + sizeof(double)*nmod; 
  mv = egl->payload + sizeof(double)* ((lmax_template+1)*(nfreq_template)) + sizeof(double)*egl->nfreq*nmod + sizeof(double)*(nmod*2);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(p=1;p<nmod;p++) {
     sprintf(name,"beamnl_I_%d_%d",(int)egl->freqlist[m1],p);
     Ip[m1*nmod+p] =  parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,);
    }
  }

  lmin =  egl->lmin;
  lmax =  egl->lmax;
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      for(ell=lmin;ell<=lmax;ell++) {
        Rq[IDX_R(egl,ell,m1,m2)] = 1;
        t1 = template[ell+(lmax_template+1)*mv[m1]];
        t2 = template[ell+(lmax_template+1)*mv[m2]];
        lt1 = log(t1);
        lt2 = log(t2);
        tp1[0] = 1./t1;
        tp2[0] = 1./t2;

        for(p=1;p<nmod;p++) {
          tp1[p] = tp1[p-1] * lt1/p;
          tp2[p] = tp2[p-1] * lt2/p;
          //_DEBUGHERE_("%d %d %d %d %g %g %g %g,%g, %g",ell,m1,m2,p,lt1,lt2,t1,t2,tp1[p]*(Ip[m1*nmod+p]-1),tp2[p]*(Ip[m2*nmod+p]-1));
          Rq[IDX_R(egl,ell,m1,m2)] += -tp1[p]*(Ip[m1*nmod+p]-1);  
          Rq[IDX_R(egl,ell,m1,m2)] += -tp2[p]*(Ip[m2*nmod+p]-1);  
        }
        for(p1=1;p1<nmod;p1++) {
          for(p2=1;p2<(nmod*(1-strict)+(strict)*(nmod-p1));p2++) {
            //_DEBUGHERE_("%d %d %d %d %d %g",ell,m1,m2,p1,p2,tp1[p1]*(Ip[m1*nmod+p1]-1)*tp2[p2]*(Ip[m2*nmod+p2]-1));
            Rq[IDX_R(egl,ell,m1,m2)] += tp1[p1]*(Ip[m1*nmod+p1]-1)*tp2[p2]*(Ip[m2*nmod+p2]-1);
          }
        }
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(beamnl,beamnl_init);
