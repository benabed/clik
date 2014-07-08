#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void dip_compute(parametric *egl, double *rq, error **err);

parametric *dip_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &dip_compute;
  
  parametric_set_default(egl,"dip_l",1800,err); 
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"dip_delta_l",300,err); 
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"A_dip_217",10,err); 
  forwardError(*err,__LINE__,NULL);
  
  
  return egl;
}

void dip_compute(parametric *egl, double *Rq, error **err) {
  double l_pivot,delta_l,v;
  int m1,m2,ell;
  double dip_nrm;
  int lmin,lmax;

  l_pivot = parametric_get_value(egl,"dip_l",err);
  forwardError(*err,__LINE__,);
  delta_l = parametric_get_value(egl,"dip_delta_l",err);
  forwardError(*err,__LINE__,);
  dip_nrm = parametric_get_value(egl,"A_dip_217",err);
  forwardError(*err,__LINE__,);
  lmin = l_pivot-delta_l;
  if (lmin<egl->lmin) {
    lmin =  egl->lmin;
  }
  lmax = l_pivot+delta_l;
  if (lmax>egl->lmax) {
    lmax =  egl->lmax;
  }

  for(m1=0;m1<egl->nfreq;m1++) {
    if (fabs(egl->freqlist[m1]-217)>1e-6) {
      continue;
    }
    for(ell=lmin;ell<=lmax;ell++) {
      Rq[IDX_R(egl,ell,m1,m1)] = -dip_nrm/ell/(ell+1.)*2.*M_PI * cos(((ell-l_pivot)/2./delta_l-1.)*M_PI) * cos(((ell-l_pivot)/2./delta_l-1.)*M_PI); 
    }
  }
}
CREATE_PARAMETRIC_FILE_INIT(dip,dip_init);
