#include "basic.h"

void powerlaw_triangle_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_triangle_fill_derivative(egl,iv,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_triangle_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_triangle_index",err);
  forwardError(*err,__LINE__,);

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

void powerlaw_tanh_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tanh_fill_derivative(egl,iv,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tanh_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tanh_index",err);
  forwardError(*err,__LINE__,);

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}
