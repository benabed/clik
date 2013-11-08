#include "cibsz.h"

// SZ alone
void sz_compute(parametric* egl, double *Rq, error **err);

parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  double fac;
  int l;

  // make sure l(l+1)/(2pi)*C_l template is normalized to 1 at l=3000 
  fac = template[3000-lmin_sz_template];
  for (l=lmin_sz_template;l<=lmax_sz_template;l++) {
    template[l-lmin_sz_template] /= fac;
  }
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  // Declare payload, allocate it and fill it

  egl->payload = malloc_err(sizeof(double)*(lmax_sz_template-lmin_sz_template+1 + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,template,sizeof(double)*(lmax_sz_template-lmin_sz_template+1));
  
  egl->eg_compute = &sz_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  parametric_set_default(egl,"sz_norm",4.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"sz_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}


void sz_compute(parametric* egl, double *Rq, error **err) {

  int ell,m1,m2,mell,nfreq;
  double *cl, *fnu, *A;
  double sz_norm, dell;
  
  nfreq = egl->nfreq;
  A = (double*) egl->payload;
  cl = A;
  fnu = &(A[lmax_sz_template-lmin_sz_template+1]);
  
  sz_norm = parametric_get_value(egl,"sz_norm",err);
  forwardError(*err,__LINE__,);

  // Compute SZ spectrum
  for (m1=0;m1<nfreq;m1++) {
    fnu[m1] = sz_spectrum((double)egl->freqlist[m1],PRM_NU0);
    printf("%g %g\n",egl->freqlist[m1],fnu[m1]);
  }

  // Create covariance matrix
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    dell = (double)ell;
    mell = (ell-lmin_sz_template);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
  Rq[IDX_R(egl,ell,m1,m2)] = sz_norm * 2.0*M_PI/(dell*(dell+1.0)) * cl[mell] * fnu[m1] * fnu[m2];
  Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}

// kSZ alone

void ksz_compute(parametric* egl, double *Rq, error **err);

parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int lmin_ksz_template= 0; //CHECK
  int lmax_ksz_template= 5000; // CHECK
  double fac;
  int l;

  // make sure l(l+1)/(2pi)*C_l template is normalized to 1 at l=3000 
  fac = template[3000-lmin_ksz_template];
  for (l=lmin_ksz_template;l<=lmax_ksz_template;l++) {
    template[l-lmin_ksz_template] /= fac;
  }
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  // Declare payload, allocate it and fill it

  egl->payload = malloc_err(sizeof(double)*(lmax_ksz_template-lmin_ksz_template+1),err);
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,template,sizeof(double)*(lmax_ksz_template-lmin_ksz_template+1));
  
  egl->eg_compute = &ksz_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  parametric_set_default(egl,"ksz_norm",0.0,err); // PICK YOUR FAVORITE
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ksz_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}


void ksz_compute(parametric* egl, double *Rq, error **err) {

  int ell,m1,m2,mell,nfreq;
  int lmin_ksz_template = 0; // CHECK
  int lmax_ksz_template = 5000; // CHECK
  double *cl,*A;
  double ksz_norm, dell;
  
  nfreq = egl->nfreq;
  A = (double*) egl->payload;
  cl = A;
  
  ksz_norm = parametric_get_value(egl,"ksz_norm",err);
  forwardError(*err,__LINE__,);

  // kSZ spectrum is like CMB

  // Create covariance matrix
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    dell = (double)ell;
    mell = (ell-lmin_ksz_template);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
  Rq[IDX_R(egl,ell,m1,m2)] = ksz_norm * 2.0*M_PI/(dell*(dell+1.0)) * cl[mell];
  Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}

//cib alone on the basis of the pip cib modelling
void ncib_compute(parametric *egl, double *rq, error **err);

parametric *ncib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv;
  int dreq[4];

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;


  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &ncib_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)* (3001*4*4) + sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* (3001*4*4));

  mv = egl->payload + sizeof(double)* (3001*4*4);

  for(m=0;m<egl->nfreq;m++) {
    double f;
    f = egl->freqlist[m];
    mv[m]=-1;
    for(i=0;i<4;i++) {
      //_DEBUGHERE_("%g %d",f,dreq[i]);
      if (fabs(f-dreq[i])<1e-6) {
        mv[m]=i;
        break;
      }  
    }
    testErrorRet(mv[m]==-1,-431432,"argl",*err,__LINE__,NULL);
  }    
  
  parametric_set_default(egl,"A_cib_217",70,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cib_index",-1.3,err); 
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ncib_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v;
  int m1,m2,ell;
  double nrm;
  int *mv;

  mv = egl->payload + sizeof(double)* (3001*4*4);

  template = egl->payload;
  l_pivot = 3000;

  index = parametric_get_value(egl,"cib_index",err);
  forwardError(*err,__LINE__,);

  nrm = parametric_get_value(egl,"A_cib_217",err);
  forwardError(*err,__LINE__,);


  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = nrm*pow((double) ell/l_pivot,(double) index-(-1.3));
    //_DEBUGHERE_("%d %g",ell,v);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*16+mv[m1]*4+mv[m2]]);
        Rq[IDX_R(egl,ell,m1,m2)] = template[ell*16+mv[m1]*4+mv[m2]] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}


// CIB with free amplitudes and power law spectrum

void cib_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,nrm;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"cib_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;
  nrm = l_pivot*(l_pivot+1)/2/M_PI;

  index = parametric_get_value(egl,"cib_index",err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }


  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }

  return;
}



parametric *cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &cib_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cib_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cib_index",-1.3,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"cib_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  parametric_add_derivative_function(egl,"any",&cib_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}

// CIB with free amplitudes and cross correlation parameter and power law spectrum

void cibr_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,nrm;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"cib_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;
  nrm = l_pivot*(l_pivot+1)/2/M_PI;

  index = parametric_get_value(egl,"cib_index",err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m1]);
    v = 1;
    v = parametric_get_value(egl,name,err);
    forwardError(*err,__LINE__,);
    A[m1*nfreq+m1] = v/nrm;
  }
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1+1;m2<nfreq;m2++) {
      sprintf(name,"cib_r_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v*sqrt(A[m1*nfreq+m1]*A[m2*nfreq+m2]);
      A[m2*nfreq+m1] = A[m1*nfreq+m2];
    }
  }


  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }

  return;
}

parametric *cibr_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &cibr_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cib_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cib_index",-1.3,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"cib_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m1]);
    parametric_set_default(egl,name,1,err);
    forwardError(*err,__LINE__,NULL);
    for(m2=m1+1;m2<egl->nfreq;m2++) {
      sprintf(name,"cib_r_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  return egl;
}

// cibsz ala gpe

void sz_cib_compute(parametric *egl, double *Rq, error **err) ;

    
parametric* sz_cib_common_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  egl =  parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  parametric_sz_cib_payload_init(egl,template,err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &sz_cib_compute;
  egl->eg_free = &parametric_sz_cib_payload_free;
  
  return egl;  
}  

void sz_cib_sz_init(parametric *egl,int der,error **err) {
  parametric_set_default(egl,"A_sz",4.0,err);
  forwardError(*err,__LINE__,);
  if (der ==1) {
    parametric_add_derivative_function(egl,"A_sz",&sz_cib_A_sz_derivative,err);
    forwardError(*err,__LINE__,);  
  }
}

void sz_cib_cib_init(parametric *egl,int der,error **err) {
 
  parametric_set_default(egl,"A_cib_100",0.,err); // change value
  forwardError(*err,__LINE__,);
  parametric_set_default(egl,"A_cib_217",70.0,err); // change value
  forwardError(*err,__LINE__,);
  parametric_set_default(egl,"A_cib_143",6.0,err); // change value
  forwardError(*err,__LINE__,);
  if (der ==1) {
    parametric_add_derivative_function(egl,"A_cib_100",&sz_cib_A_cib_100_derivative,err);
    forwardError(*err,__LINE__,);
    parametric_add_derivative_function(egl,"A_cib_143",&sz_cib_A_cib_143_derivative,err);
    forwardError(*err,__LINE__,);
    parametric_add_derivative_function(egl,"A_cib_217",&sz_cib_A_cib_217_derivative,err);
    forwardError(*err,__LINE__,);
  }

  parametric_set_default(egl,"r_cib_143_217",0.85,err); // change value
  forwardError(*err,__LINE__,);
  parametric_set_default(egl,"r_cib_100_143",0.,err); // change value
  forwardError(*err,__LINE__,);
  parametric_set_default(egl,"r_cib_100_217",0.,err); // change value
  forwardError(*err,__LINE__,);
  
  if (der ==1) {
    parametric_add_derivative_function(egl,"r_cib_100_143",&sz_cib_r_cib_100_143_derivative,err);
    forwardError(*err,__LINE__,);
    parametric_add_derivative_function(egl,"r_cib_100_217",&sz_cib_r_cib_100_217_derivative,err);
    forwardError(*err,__LINE__,);
    parametric_add_derivative_function(egl,"r_cib_143_217",&sz_cib_r_cib_143_217_derivative,err);
    forwardError(*err,__LINE__,);
  }

  parametric_set_default(egl,"cib_index",-1.3,err); // change value
  forwardError(*err,__LINE__,);
  if (der ==1) {
    parametric_add_derivative_function(egl,"cib_index",&sz_cib_cib_index_derivative,err);
    forwardError(*err,__LINE__,);
  }
}

void sz_cib_szxcib_init(parametric *egl,int der,error **err) {
  sz_cib_payload *payload;
  long remove_100;
  int m1,m2;

  payload = egl->payload;

  remove_100 = 0;
  remove_100 = pflist_get_int_value(egl->pf,"no_szxcib_100",&remove_100,err);
  forwardError(*err,__LINE__,);

  if (remove_100==1) {
    payload = egl->payload;
    for(m1=0;m1<egl->nfreq;m1++) {
      if(payload->ind_freq[m1]==0) {
        payload->ind_freq[m1]=-1;
      }
    }
  }
  parametric_set_default(egl,"xi_sz_cib",0.5,err); // change value
  forwardError(*err,__LINE__,);
  if (der ==1) {
    parametric_add_derivative_function(egl,"xi_sz_cib",&sz_cib_xi_sz_cib_derivative,err);
    forwardError(*err,__LINE__,);
  }
}

parametric *sz_cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  
  // Initialize structure and payload, with template
  egl = sz_cib_common_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,template,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_sz_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_cib_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_szxcib_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);
  
  pflist_set_int_value(egl->pf,"do_szxcib",1,err);
  pflist_set_int_value(egl->pf,"do_cib",1,err);
  pflist_set_int_value(egl->pf,"do_sz",1,err);
    
  return egl;

}

parametric *sz_cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  
  // Initialize structure and payload, with template
  egl = sz_cib_common_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,template,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_sz_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_cib_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_szxcib_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);
  
  pflist_set_int_value(egl->pf,"do_szxcib",1,err);
  pflist_set_int_value(egl->pf,"do_cib",0,err);
  pflist_set_int_value(egl->pf,"do_sz",0,err);
    
  return egl;

}

parametric *cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax,  error **err) {
  parametric *egl;
  
  // Initialize structure and payload, with template
  egl = sz_cib_common_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,NULL,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_cib_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);
  
  pflist_set_int_value(egl->pf,"do_szxcib",0,err);
  pflist_set_int_value(egl->pf,"do_cib",1,err);
  pflist_set_int_value(egl->pf,"do_sz",0,err);
    
  return egl;
}

parametric *sz_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  
  // Initialize structure and payload, with template
  egl = sz_cib_common_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,template,err);
  forwardError(*err,__LINE__,NULL);

  sz_cib_sz_init(egl,1,err);
  forwardError(*err,__LINE__,NULL);

  pflist_set_int_value(egl->pf,"do_szxcib",0,err);
  pflist_set_int_value(egl->pf,"do_cib",0,err);
  pflist_set_int_value(egl->pf,"do_sz",1,err);
    
  return egl;
}


void sz_cib_compute(parametric *egl, double *Rq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
  
  parametric_zero_rq(egl,Rq);

  // Compute the SZ part first
  if (do_sz ==1) { 

    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = (ell-lmin_sz_template);
      for (m1=0;m1<egl->nfreq;m1++) {
        for (m2=m1;m2<egl->nfreq;m2++) {
          Rq[IDX_R(egl,ell,m1,m2)] += a_sz * 2.0*M_PI/(dell*(dell+1.0)) * sz_template[mell] * fnu[m1] * fnu[m2];
          Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
        }
      }
    }
  }
  
  // Add the CIB part
  if (do_cib ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {      
      for(m1=0;m1<nfreq;m1++) {
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++){
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            Rq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind1]*a_cib[ind2])/d3000 * pow(((double)ell/(double)lnorm),cib_index);
            Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
          }        
        }
      }
    }
  }
  
  // Add the CIB-SZ correlation part
  if (do_szxcib ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            Rq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * ( sqrt(fnu[m1]*a_cib[ind2]) + sqrt(fnu[m2]*a_cib[ind1]) ) *
               corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
            Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }
    
  return;
}  

CREATE_PARAMETRIC_FILE_INIT(cib,cib_init);
CREATE_PARAMETRIC_FILE_INIT(cibr,cibr_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz,sz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ksz,ksz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_cib,sz_cib_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_x,sz_x_init);
CREATE_PARAMETRIC_FILE_INIT(cib_x,cib_x_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_cib_x,sz_cib_x_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ncib,ncib_init);

