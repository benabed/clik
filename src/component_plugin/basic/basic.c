#include "clik_parametric.h"
#include "clik_parametric_addon.h"

#define PRM_NU0 143.

// some helpfunctions

void basic_compute_AB(parametric *egl, error **err) {
  int m1,m2,nfreq;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  A = egl->payload;

  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  for (m1=0;m1<nfreq;m1++) {
    vec[m1] = dBdT((double)egl->freqlist[m1],PRM_NU0);
    for(m2=m1;m2<nfreq;m2++) {
      x = (double)egl->freqlist[m1]*(double)egl->freqlist[m2]/(PRM_NU0*PRM_NU0);
      lnx = log(x);
      ln2x = lnx*lnx;
      A[m1*nfreq+m2] = lnx;
      A[m2*nfreq+m1] = lnx;
      B[m1*nfreq+m2] = ln2x;
      B[m2*nfreq+m1] = ln2x;
    }
  }
  return;
}


void basic_alpha_derivative(parametric *egl,int iv, double *Rq, double *dRq, error **err) {
  double *A;
  int ell,m1,m2;

  A = egl->payload;
  // dR/dalpha_rg = log(nu1*nu2/nu0^2) * R
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * Rq[IDX_R(egl,ell,m1,m2)];
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}

void basic_sigma_derivative(parametric *egl,int iv, double *Rq, double *dRq, error **err) {
  double *A,*B;
  double sigma;
  int ell,m1,m2;
  
  A = egl->payload;
  B = A + egl->nfreq*egl->nfreq;
  
  sigma = parametric_get_value(egl,egl->varkey[iv],err);
  forwardError(*err,__LINE__,);
  // dR/dsigma_rg = log(nu1*nu2/nu0^2)^2 * R
  
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = sigma * B[m1*egl->nfreq+m2] * Rq[IDX_R(egl,ell,m1,m2)];
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}

// Millea et al. radio galaxies Poisson contribution. Based on a single population with given mean spectral index, as well as dispersion

void radiogal_compute(parametric* exg, double *Rq, error **err);

parametric *radiogal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &radiogal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"radiogal_norm",78.5,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"radiogal_alpha",-0.36,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"radiogal_sigma",0.64,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "radiogal_sigma", &basic_sigma_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void radiogal_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm_rg, alpha_rg, sigma_rg;
  double d3000;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  d3000 = 3000.*3001./2./M_PI;
  
  A = egl->payload;
  
  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  
  norm_rg = parametric_get_value(egl,"radiogal_norm",err);
  forwardError(*err,__LINE__,);

  alpha_rg = parametric_get_value(egl,"radiogal_alpha",err);
  forwardError(*err,__LINE__,);

  sigma_rg = parametric_get_value(egl,"radiogal_sigma",err);
  forwardError(*err,__LINE__,);

  basic_compute_AB(egl, err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = norm_rg/d3000 * exp(alpha_rg*A[m1*nfreq+m2] +
              sigma_rg*sigma_rg/2.0 * B[m1*nfreq+m2])/(vec[m1]*vec[m2]);
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////
// CIB model. Millea et al. for now
////////////////////////////////////////////////////////////////////////////////////////////

void ir_poisson_compute(parametric* exg, double *Rq, error **err);

parametric *ir_poisson_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &ir_poisson_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_norm",5.9,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_alpha",3.8,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_poisson_sigma",0.4,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_poisson_sigma", &basic_sigma_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void ir_poisson_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_poisson_norm, ir_poisson_alpha, ir_poisson_sigma;
  double d3000;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  basic_compute_AB(egl,err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  d3000 = 3000.*3001./2./M_PI;

  ir_poisson_norm = parametric_get_value(egl,"ir_poisson_norm",err);
  forwardError(*err,__LINE__,);

  ir_poisson_alpha = parametric_get_value(egl,"ir_poisson_alpha",err);
  forwardError(*err,__LINE__,);

  ir_poisson_sigma = parametric_get_value(egl,"ir_poisson_sigma",err);
  forwardError(*err,__LINE__,);

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = ir_poisson_norm/d3000 * exp(ir_poisson_alpha*A[m1*nfreq+m2] +
              ir_poisson_sigma*ir_poisson_sigma/2.0 * B[m1*nfreq+m2])/(vec[m1]*vec[m2]);
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////
// CIB clustered part. Millea et al. for now.
////////////////////////////////////////////////////////////////////////////////////////

void ir_clustered_compute(parametric* exg, double *Rq, error **err);
void ir_clustered_step_derivative(parametric* egl, int iv, double *Rq, double *dRq, error **err);

parametric *ir_clustered_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  pfchar type;
  char *pt;


  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &ir_clustered_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq + (lmax-lmin+1)),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_norm",3.9,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_norm", &parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_alpha",3.8,err); // Millea et al. ref value
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_alpha", &basic_alpha_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ir_clustered_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "ir_clustered_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(type,"step");
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"ir_clustered_correlation",pt,err);
  forwardError(*err,__LINE__,NULL);

  //_DEBUGHERE_("%s",pt);
  
  if (strcmp(pt,"step")==0) {
    parametric_set_default(egl,"ir_clustered_correlation_step",0.5,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl, "ir_clustered_correlation_step", &ir_clustered_step_derivative,err);
    forwardError(*err,__LINE__,NULL);
  } else if (strcmp(pt,"matrix")==0) {
    int m1,m2;
    double vl;
    pfchar name;
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        sprintf(name,"ir_clustered_correlation_M_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
        parametric_set_default(egl,name,1,err);
        forwardError(*err,__LINE__,NULL);
      }
    }  
  } else {
    testErrorRetVA(1==1,-1234,"Unknown ir_clustered_correlation type '%s'",*err,__LINE__,,type);
    // return ?
  }

  return egl;
}

double ir_clustered_step_index(parametric* egl, int m1,int m2) {
  return fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33;
}

void ir_clustered_step_derivative(parametric* egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2; 
  double step ;
  char *key;

  key = egl->varkey[iv];
  step = parametric_get_value(egl,key,err);
  forwardError(*err,__LINE__,);
  
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      dRq[IDX_R(egl,ell,m1,m1)] = 0; //diagonal unchanged
      for (m2=m1+1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = Rq[IDX_R(egl,ell,m1,m2)]/step * ir_clustered_step_index(egl,m1,m2);
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}

void ir_clustered_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double ir_clustered_norm, ir_clustered_alpha,ir_clustered_index;
  double d3000,t3000;
  double x, lnx, ln2x;
  double *A,*vec,*template,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;

  A = egl->payload; // Will store log(nu/nu0)
  d3000 = 3000.*3001./2./M_PI;
  
  egl->l_pivot = 3000;

  ir_clustered_norm = parametric_get_value(egl,"ir_clustered_norm",err);
  forwardError(*err,__LINE__,);

  ir_clustered_alpha = parametric_get_value(egl,"ir_clustered_alpha",err);
  forwardError(*err,__LINE__,);

  ir_clustered_index = parametric_get_value(egl,"ir_clustered_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  vec = A + nfreq*nfreq; // Will store dB/dT(nu,nu0)
  template = vec + nfreq; // Will store Cl spectrum (unnormalized)

  for (m1=0;m1<nfreq;m1++) {
    vec[m1] = dBdT((double)egl->freqlist[m1],PRM_NU0);
    for(m2=m1;m2<nfreq;m2++) {
      x = (double)egl->freqlist[m1]*(double)egl->freqlist[m2]/(PRM_NU0*PRM_NU0);
      lnx = log(x);
      A[m1*nfreq+m2] = lnx;
      A[m2*nfreq+m1] = lnx;
    }
  }

  // Create Cl template
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=ell-egl->lmin;
    template[mell] = pow((double)ell,ir_clustered_index);
  }
  t3000 = pow(3000.,ir_clustered_index);

  //create correlation matrix
  dcm = template + egl->lmax-egl->lmin+1;

  sprintf(type,"step");
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"ir_clustered_correlation",pt,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("%s",pt);

  isstep=0;
  if (strcmp(pt,"step")==0) {  
    isstep=1;

    step = parametric_get_value(egl,"ir_clustered_correlation_step",err);
    forwardError(*err,__LINE__,);
    for(m1=0;m1<egl->nfreq;m1++) {
      dcm[m1*egl->nfreq+m1] = 1.;
      for(m2=0;m2<m1;m2++) {
        dcm[m1*egl->nfreq+m2] = pow(step,ir_clustered_step_index(egl,m1,m2));
        dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
        //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
      }
    }
  } else {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        pfchar name;
        sprintf(name,"ir_clustered_correlation_M_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
        dcm[m1*egl->nfreq+m2] = parametric_get_value(egl,name,err);
        forwardError(*err,__LINE__,);
        dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
        //_DEBUGHERE_("%d %d %g",m1,m2,dcm[m1*egl->nfreq+m2]);
      }
    }
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    lell = ell - egl->lmin;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = ir_clustered_norm * 
          template[lell]/(d3000*t3000) * 
          exp(ir_clustered_alpha*A[m1*nfreq+m2])/(vec[m1]*vec[m2]) * dcm[m1*egl->nfreq+m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}




////////////////////////////////////////////////////////////////////////////////////////
// Galactic model
////////////////////////////////////////////////////////////////////////////////////////

// Dust emissivity normalized to 1 at nu0
double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res, res0;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  res0 = pow(nu0,3.0+beta_dust)/(ex0 - 1.); // Overall normalization will go away
  x = nu * h_over_kT;
  ex = exp(x);
  res = pow(nu,3.0+beta_dust)/(ex - 1.);

  // Return dust emissivity normalized to one at nu0, in dT (CMB)
  return (res/res0/dBdT(nu,nu0));
}

// Derivative of dust_spectrum with respect to T_dust
double d_dust_spectrum_d_T_dust(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  x = nu * h_over_kT;
  ex = exp(x);
  res = 1./T_dust * pow(nu/nu0,3.0+beta_dust) * (x*ex*(ex0-1.)-x0*ex0*(ex-1.)) / ((ex-1.)*(ex-1));
  res /= dBdT(nu,nu0);

  return(res);
  
} 

// Derivative of dust_spectrum with respect to beta_dust
double d_dust_spectrum_d_beta_dust(double nu, double T_dust, double beta_dust, double nu0) {
  
  return (log(nu/nu0)*dust_spectrum(nu,T_dust,beta_dust,nu0));

}

// Non-thermal emissivity normalized to 1 at nu0 (e.g. synchrotron, free-free, etc.)
// The spectral index is called alpha here, to distinguish from the gray body case
// Expressed in dT (CMB)
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0) {

  //return (pow(nu/nu0,alpha_non_thermal)/dBdT(nu,nu0));
  return (exp(alpha_non_thermal*log(nu/nu0))/dBdT(nu,nu0));
}

// Derivative of non_thermal_spectrum with respect to alpha
double d_non_thermal_spectrum_d_alpha_non_thermal(double nu, double alpha_non_thermal, double nu0) {

  return (log(nu/nu0)*non_thermal_spectrum(nu,alpha_non_thermal,nu0));
}

void galactic_component_compute(parametric* exg, double *Rq, error **err);


void galactic_component_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm,l_pivot,index,alpha_non_thermal,beta_dust,T_dust;
  double v,lA;
  double *A, *a, *B, *b;
  pfchar type;
  char* pt;
  int isdust;
  double prm_nu0;
  
  // dust or non_thermal

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  if (strcmp(pt,"dust")==0) {
    
    beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
    forwardError(*err,__LINE__,);

    T_dust = parametric_get_value(egl,"gal_T_dust",err);
    forwardError(*err,__LINE__,);

  } else if (strcmp(pt,"non_thermal")==0) {

    isdust = 0;
    alpha_non_thermal = parametric_get_value(egl,"gal_alpha_non_thermal",err);
    forwardError(*err,__LINE__,);

  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,,type);
  }
    

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  for (m1=0;m1<nfreq;m1++) {
    if (isdust) {
      a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
    } else {
      a[m1] = non_thermal_spectrum((double)egl->freqlist[m1],alpha_non_thermal,prm_nu0);
    }
  }

  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      v = a[m1]*a[m2];
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  // R_q

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index) * norm;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }

  return;
}
void gal_beta_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double T_dust, beta_dust;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_T_dust",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dbeta_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_dust_spectrum_d_beta_dust((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;        
}

void gal_T_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double T_dust, beta_dust;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_T_dust",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dT_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_dust_spectrum_d_T_dust((double)egl->freqlist[m1],T_dust,beta_dust,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
}

void gal_alpha_non_thermal_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  int ell,m1,m2;
  double *A, *B, *a, *b;
  int nfreq;
  double alpha_non_thermal;
  double lA,v;
  double norm,l_pivot,index;
  double prm_nu0;

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  alpha_non_thermal = parametric_get_value(egl,"gal_alpha_non_thermal",err);
  forwardError(*err,__LINE__,);
  
  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);
  
  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  prm_nu0 = parametric_get_value(egl,"gal_nu0",err);
  forwardError(*err,__LINE__,);

  // dR/dalpha_non_thermal
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    b[m1] = d_non_thermal_spectrum_d_alpha_non_thermal((double)egl->freqlist[m1],alpha_non_thermal,prm_nu0);
  }
  for (m1=0;m1<nfreq;m1++) {
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = b[m1]*a[m2] + b[m2]*a[m1];
      B[m2*nfreq+m1] = B[m1*nfreq+m2];
    }
  }
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/egl->l_pivot,index) * norm;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        lA = B[m1*nfreq+m2];
        dRq[IDX_R(egl,ell,m1,m2)] = lA*v;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }  
}

parametric *galactic_component_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  pfchar type;
  char* pt;
  int isdust;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &galactic_component_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*2*egl->nfreq*(egl->nfreq+1),err);
  forwardError(*err,__LINE__,NULL);
  
  // uK^2 at l=500, nu=143 GHz;
  parametric_set_default(egl,"gal_norm",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_norm",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_norm",&parametric_norm_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_index",0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_index",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,NULL);

  if (strcmp(pt,"dust")==0) {
    parametric_set_default(egl,"gal_beta_dust",1.8,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_beta_dust",&gal_beta_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

    parametric_set_default(egl,"gal_T_dust",18,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_T_dust",&gal_T_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

  } else if (strcmp(pt,"non_thermal")==0) {
    //Intensity, = -3.0 in RJ
    parametric_set_default(egl,"gal_alpha_non_thermal",-1,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_alpha_non_thermal",&gal_alpha_non_thermal_derivative,err);  
    forwardError(*err,__LINE__,NULL);
  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,,type);
    // return ?
  }

  return egl;
}

void poisson_tensor_bydet_compute(parametric* exg, double *Rq, error **err);

parametric *poisson_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  
  egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &poisson_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","poisson_tensor_bydet_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    double nrm,frq;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    frq = egl->freqlist[ic];
    nrm = sqrt(80/(3000.*3001./2./M_PI)*(1.98984227 * frq*frq*frq -27.70516974 *frq*frq + 125.04493152 *frq -181.36576083));
    parametric_set_default(egl,Ac,nrm,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;

  return egl;
}

void poisson_tensor_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double *A;
  int ic;
  pfchar Ac;

  A = egl->payload;

  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  
  
  nfreq = egl->nfreq;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  
  return;
}

void powerlaw_tensor_bydet_compute(parametric* exg, double *Rq, error **err);

parametric *powerlaw_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_bydet_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_bydet_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_bydet_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_bydet_correlation_step",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_bydet_correlation_step", &ir_clustered_step_derivative,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","powerlaw_tensor_bydet_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  parametric_set_default(egl,"powerlaw_tensor_bydet_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;
  return egl;
}

void powerlaw_tensor_bydet_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tensor_bydet_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tensor_bydet_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  step = parametric_get_value(egl,"powerlaw_tensor_bydet_correlation_step",err);  
  forwardError(*err,__LINE__,);

  dcm = A + egl->nfreq;

  for(m1=0;m1<egl->nfreq;m1++) {
    dcm[m1*egl->nfreq+m1] = 1.;
    for(m2=0;m2<m1;m2++) {
      dcm[m1*egl->nfreq+m2] = pow(step,ir_clustered_step_index(egl,m1,m2));
      dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
      //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1]*A[m2] * dcm[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

parametric *poisson_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;
  
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &poisson_tensor_bydet_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*(egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","poisson_tensor_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    double nrm,frq;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    frq = egl->freqlist[ic];
    nrm = sqrt(80/(3000.*3001./2./M_PI)*(1.98984227 * frq*frq*frq -27.70516974 *frq*frq + 125.04493152 *frq -181.36576083));
    parametric_set_default(egl,Ac,nrm,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;

  return egl;
}



void powerlaw_tensor_compute(parametric* exg, double *Rq, error **err);

parametric *powerlaw_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tensor_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*(egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tensor_correlation_step",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tensor_correlation_step", &ir_clustered_step_derivative,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(egl->tensor_norm_template,"%s","powerlaw_tensor_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  parametric_set_default(egl,"powerlaw_tensor_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &parametric_tensor_norm_derivative;
  return egl;
}

void powerlaw_tensor_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tensor_fill(egl,A,err);
  forwardError(*err,__LINE__,);  

  l_pivot = parametric_get_value(egl,"powerlaw_tensor_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"powerlaw_tensor_index",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  step = parametric_get_value(egl,"powerlaw_tensor_correlation_step",err);  
  forwardError(*err,__LINE__,);

  dcm = A + egl->nfreq;

  for(m1=0;m1<egl->nfreq;m1++) {
    dcm[m1*egl->nfreq+m1] = 1.;
    for(m2=0;m2<m1;m2++) {
      dcm[m1*egl->nfreq+m2] = pow(step,ir_clustered_step_index(egl,m1,m2));
      dcm[m2*egl->nfreq+m1] = dcm[m1*egl->nfreq+m2];
      //_DEBUGHERE_("%d %d %g %g",m1,m2,dcm[m1*egl->nfreq+m2],fabs(log(egl->freqlist[m1])-log(egl->freqlist[m2]))*2.33)
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1]*A[m2] * dcm[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

void powerlaw_triangle_compute(parametric* exg, double *Rq, error **err);
void powerlaw_triangle_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);


parametric *powerlaw_triangle_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_triangle_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*((egl->nfreq*egl->nfreq)*2),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_triangle_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_triangle_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  
  sprintf(egl->tensor_norm_template,"%s","powerlaw_triangle_T_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    int jc;
    for(jc=ic;jc<egl->nfreq;jc++) {
      sprintf(Ac,"%s%d_%d",egl->tensor_norm_template,ic,jc);    
      parametric_set_default(egl,Ac,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }
  
  parametric_set_default(egl,"powerlaw_triangle_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &powerlaw_triangle_norm_derivative;
  return egl;
}

void powerlaw_triangle_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_triangle_fill(egl,A,err);
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
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
  return;
}

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


void powerlaw_tanh_compute(parametric* exg, double *Rq, error **err);
void powerlaw_tanh_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);


parametric *powerlaw_tanh_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int ic;
  pfchar Ac;

  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_tanh_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)*((egl->nfreq*egl->nfreq)*2),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"powerlaw_tanh_index",-1.25,err); // Karim's favorite
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl, "powerlaw_tanh_index", &parametric_index_derivative,err);
  forwardError(*err,__LINE__,NULL);

  
  sprintf(egl->tensor_norm_template,"%s","powerlaw_tanh_A_");
  egl->tensor_norm_template_len = strlen(egl->tensor_norm_template);
  
  for(ic=0;ic<egl->nfreq;ic++) {
    int jc;
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    parametric_set_default(egl,Ac,1,err);
    forwardError(*err,__LINE__,NULL);
    for(jc=ic+1;jc<egl->nfreq;jc++) {
      sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,ic,jc);
      parametric_set_default(egl,Ac,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }
  
  parametric_set_default(egl,"powerlaw_tanh_l_pivot",1000,err);
  forwardError(*err,__LINE__,NULL);
    
  egl->eg_deriv_any = &powerlaw_tanh_norm_derivative;
  return egl;
}

void powerlaw_tanh_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,lell;
  double l_pivot,index,v;
  double x, lnx, ln2x;
  double *A,*dcm;
  pfchar type;
  char *pt;
  double step;
  int isstep;
  
  A = egl->payload;
  parametric_tanh_fill(egl,A,err);
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
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * v;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
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


// SZ alone
void sz_compute(parametric* egl, double *Rq, error **err);

parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int lmin_sz_template= 2;
  int lmax_sz_template= 10000;
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
  int lmax_sz_template = 10000;
  int lmin_sz_template = 2;
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

// SZ+CIB+CROSS

#define N_FREQ_CIB 3
int CIB_FREQ[] = {100,143,217};
double CIB_FIX_COLOR[] = {1.122, 1.134, 1.33};
double SQ_SZ_FIX_COLOR[] = {0.9619, 0.95, 0.0};
#define lmin_sz_template 2
#define lmax_sz_template  10000                               
#define lmin_corr_template  2
#define lmax_corr_template  9999

typedef struct {
  int *ind_freq;
  double *template;
  double *sz_template;
  double *corr_template;
  double a_cib[N_FREQ_CIB];
  double r_cib[N_FREQ_CIB*N_FREQ_CIB];
  double *fnu;
} sz_cib_payload;

void parametric_sz_cib_payload_free(void **pp) {
  sz_cib_payload *p;
  p = *pp;
  if (p!=NULL) {
    free(p->template);
    free(p->ind_freq);
    free(p);
  }
  *pp=NULL;
}


void parametric_sz_cib_payload_init(parametric *egl, double *template, error **err) {
  int m1,m2;
  int l;
  double fac;
  int lnorm = 3000;
  sz_cib_payload *payload;
  int template_size;
  long fix_sz_color;

  testErrorRetVA(egl->lmax>lmax_sz_template,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,,egl->lmax,lmax_sz_template);
  testErrorRetVA(egl->lmax>lmax_corr_template,-1111,"lmax too big (got %d should be < %d",*err,__LINE__,,egl->lmax,lmax_corr_template);

  egl->payload = malloc_err(sizeof(sz_cib_payload), err);
  forwardError(*err, __LINE__,);

  payload = egl->payload;
  
  payload ->ind_freq = malloc_err(sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,);
  
  for (m1=0;m1<egl->nfreq;m1++) {
    payload->ind_freq[m1]=-1;
    for (m2=0;m2<N_FREQ_CIB;m2++) {
      if (fabs(egl->freqlist[m1]-CIB_FREQ[m2])<1e-6) {
         payload->ind_freq[m1]=m2;
      }
    }
  }

  template_size = lmax_sz_template - lmin_sz_template + 1 + lmax_corr_template - lmin_corr_template + 1;
  payload->template = malloc_err(sizeof(double)*(template_size+egl->nfreq),err);
  forwardError(*err,__LINE__,);
  if (template!=NULL) {
    memcpy(payload->template,template,sizeof(double)*template_size);
    payload->sz_template = payload->template;
    payload->corr_template = payload->template + lmax_sz_template - lmin_sz_template + 1;
    fac = payload->sz_template[lnorm-lmin_sz_template];
    for (l=lmin_sz_template;l<=lmax_sz_template;l++) {
      payload->sz_template[l-lmin_sz_template] /= fac;
    }
    fac = payload->corr_template[lnorm-lmin_corr_template];
    for (l=lmin_corr_template;l<=lmax_corr_template;l++) {
      payload->corr_template[l-lmin_corr_template] /= fac;
    }
  }
  payload->fnu = payload->template + template_size;
  fix_sz_color = 1;
  fix_sz_color = pflist_get_int_value(egl->pf,"use_sz_fix_color",&fix_sz_color,err);             
  forwardError(*err,__LINE__,);     
  for (m1=0;m1<egl->nfreq;m1++) {                                        
    payload->fnu[m1] = sz_spectrum((double)egl->freqlist[m1],PRM_NU0);  
    if (fix_sz_color !=0) {
      payload->fnu[m1] *= sqrt(SQ_SZ_FIX_COLOR[payload->ind_freq[m1]]);
    }
  }
  for(m1=0;m1<N_FREQ_CIB;m1++) {
    payload->r_cib[m1*N_FREQ_CIB+m1] = 1;
    for(m2=m1+1;m2<N_FREQ_CIB;m2++) {
      payload->r_cib[m1*N_FREQ_CIB+m2] = 0;
      payload->r_cib[m2*N_FREQ_CIB+m1] = 0;
    }
  }                                                                 
}

void sz_cib_compute(parametric *egl, double *Rq, error **err) ;
void sz_cib_A_sz_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_100_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_A_cib_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_100_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_100_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_r_cib_143_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_xi_sz_cib_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void sz_cib_cib_index_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);

void sz_cib_fill_cib_struct(parametric *egl,error **err) {
  int i,j;
  char tmpl[1000];
  long fix_color;
  sz_cib_payload *payload;
  
  payload = egl->payload;

  fix_color = 1;
  fix_color = pflist_get_int_value(egl->pf,"use_cib_fix_color",&fix_color,err);             
  forwardError(*err,__LINE__,);                                       
  
  
  for(i=0;i<N_FREQ_CIB;i++) {
    sprintf(tmpl,"A_cib_%d",CIB_FREQ[i]);
    payload->a_cib[i] = parametric_get_value(egl,tmpl,err);
    forwardError(*err,__LINE__,);
    if (fix_color !=0) {
      payload->a_cib[i] *= CIB_FIX_COLOR[i];
    }
    for(j=i+1;j<N_FREQ_CIB;j++) {
      sprintf(tmpl,"r_cib_%d_%d",CIB_FREQ[i],CIB_FREQ[j]);
      payload->r_cib[i*N_FREQ_CIB+j] = parametric_get_value(egl,tmpl,err);
      forwardError(*err,__LINE__,);
      payload->r_cib[j*N_FREQ_CIB+i] =payload->r_cib[i*N_FREQ_CIB+j];
    }
  }
}

#define SZ_CIB_DEFS                          \
  double a_sz,  xi_sz_cib;                   \
  int lnorm = 3000;                          \
  double *sz_template, *corr_template, *fnu; \
  double *a_cib,*r_cib;                      \
  double cib_corr, cib_index;                \
  sz_cib_payload *payload;                   \
  int do_sz,do_cib,do_szxcib;                \
  int ell,m1,m2,mell,nfreq;                  \
  int *ind_freq;                             \
  int ind1, ind2;                            \
  double d3000, dell;
  

#define SZ_CIB_INITS                                                  \
  do_sz = pflist_get_int_value(egl->pf,"do_sz",NULL,err);             \
  forwardError(*err,__LINE__,);                                       \
  do_cib = pflist_get_int_value(egl->pf,"do_cib",NULL,err);           \
  forwardError(*err,__LINE__,);                                       \
  do_szxcib = pflist_get_int_value(egl->pf,"do_szxcib",NULL,err);     \
  forwardError(*err,__LINE__,);                                       \
  nfreq = egl->nfreq;                                                 \
  payload = egl->payload;                                             \
  ind_freq = payload->ind_freq;                                       \
  d3000 = (3000.0*3001.)/2.0/M_PI;                                    \
  if (do_cib == 1 || do_szxcib == 1) {                                \
    a_cib = payload->a_cib;                                           \
    r_cib = payload->r_cib;                                           \
    sz_cib_fill_cib_struct(egl,err);                                  \
    forwardError(*err,__LINE__,);                                     \
    cib_index = parametric_get_value(egl,"cib_index",err);            \
    forwardError(*err,__LINE__,);                                     \
  }                                                                   \
  if (do_szxcib == 1)  {                                              \
    xi_sz_cib = parametric_get_value(egl,"xi_sz_cib",err);            \
    forwardError(*err,__LINE__,);                                     \
    sz_template = payload->sz_template;                               \
    corr_template = payload->corr_template;                           \
  }                                                                   \
  if (do_sz == 1 || do_szxcib == 1)  {                                \
    a_sz = parametric_get_value(egl,"A_sz",err);                      \
    forwardError(*err,__LINE__,);                                     \
    sz_template = payload->sz_template;                               \
    fnu = payload->fnu;                                               \
  }                                                                    
  

// Apply frequency averaged color corrections
// Numerical values set to agree with Camspec

  
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

void parametric_zero_rq(parametric *egl, double* rrq) {
  int ell,m1,m2;
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {  
    for(m1=0;m1<egl->nfreq;m1++) {                
      for(m2=m1;m2<egl->nfreq;m2++) {             
        rrq[IDX_R(egl,ell,m1,m2)] = 0;       
        rrq[IDX_R(egl,ell,m2,m1)] = 0;       
      }                                      
    }                                        
  }                                          
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


void sz_cib_A_sz_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
    
  parametric_zero_rq(egl,dRq);
  
  // Compute the SZ part first
  if (do_sz ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = (ell-lmin_sz_template);
      for (m1=0;m1<nfreq;m1++) {
        for (m2=m1;m2<nfreq;m2++) {
          dRq[IDX_R(egl,ell,m1,m2)] += 2.0*M_PI/(dell*(dell+1.0)) * sz_template[mell] * fnu[m1] * fnu[m2];
          dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
        }
      }
    }
  }

  // NO CIB part
  // Add the CIB-SZ correlation part
  if (do_szxcib ==1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz for now
            dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * 1./(2.0*sqrt(a_sz)) * ( sqrt(fnu[m1]*a_cib[ind2]) + sqrt(fnu[m2]*a_cib[ind1]) ) *
               corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }
  return;
}  

void sz_cib_A_cib_XXX_derivative(parametric *egl, int II,int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
  
  parametric_zero_rq(egl,dRq);
  
  // NO SZ part
  // Add the CIB part
  if (do_cib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      for(m1=0;m1<nfreq;m1++) {
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++){
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            if ((ind1==II) || (ind2==II)) { // freq1 or freq2 = 100 GHz
              if (ind1==ind2) { // freq1 = freq2 = 100 GHz
                dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]/d3000 * pow(((double)ell/(double)lnorm),cib_index);
              } else { // freq1 != freq2 = 143, 217 GHz
                if (ind1==II) {
                  dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind2]/a_cib[ind1])/2.0/d3000 * pow(((double)ell/(double)lnorm),cib_index);  
                } else {
                  dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind1]/a_cib[ind2])/2.0/d3000 * pow(((double)ell/(double)lnorm),cib_index);
                }
                dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
              }
            }
          }
        }
      }
    }
  }
  
  // Add the CIB-SZ correlation part
  if (do_szxcib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      dell = (double)ell;
      mell = ell - lmin_corr_template;
      for (m1=0;m1<nfreq;m1++){
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++) {
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            if ((ind1==II) || (ind2==II)) { // freq1 or freq2 = 100 GHz
              if (ind1==ind2) { // freq1 = freq2 = 100 GHz
                dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m1]/a_cib[ind1]) * corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
              } else { 
                if (ind1==II) {
                  dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m2]/a_cib[ind1])/2.0 *  corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
                } else {
                  dRq[IDX_R(egl,ell,m1,m2)] -= xi_sz_cib * sqrt(a_sz) * sqrt(fnu[m1]/a_cib[ind2])/2.0 *  corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
                }
                dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
              }
            }
          }
        }
      }
    }
  }
  return;
}

void sz_cib_A_cib_100_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,0,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_A_cib_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,1,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_A_cib_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_A_cib_XXX_derivative(egl,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_XXX_YYY_derivative(parametric *egl, int II, int JJ, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;
  
  parametric_zero_rq(egl,dRq);

  /* NO SZ part*/                                                                                                          
  /* Add the CIB part */                                                                                                   
  if (do_cib == 1)  {                                                                                                      
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {                                                                             
      for(m1=0;m1<nfreq;m1++) {                                                                                            
        ind1 = ind_freq[m1];
        if (ind1!=II && ind1!=JJ) {
          continue;
        }                                                                                               
        for (m2=m1;m2<nfreq;m2++){                                                                                         
          ind2 = ind_freq[m2];                                                                                             
          if (((ind1 == II) && (ind2 == JJ)) || ((ind1 == JJ) && (ind2 == II))) {                                                                              
            dRq[IDX_R(egl,ell,m1,m2)] += sqrt(a_cib[ind1]*a_cib[ind2])/d3000 * pow(((double)ell/(double)lnorm),cib_index); 
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];                                                         
          }                                                                                                                
        }                                                                                                                  
      }                                                                                                                    
    }                                                                                                                      
  }                                                                                                                        
  return; 
}

void sz_cib_r_cib_100_143_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_r_cib_XXX_YYY_derivative(egl,0,1,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_100_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  sz_cib_r_cib_XXX_YYY_derivative(egl,0,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_r_cib_143_217_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
sz_cib_r_cib_XXX_YYY_derivative(egl,1,2,iv,Rq,dRq,err);
  forwardError(*err,__LINE__,);
  return;
}

void sz_cib_xi_sz_cib_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;  

  parametric_zero_rq(egl,dRq);

  // NO SZ part, no CIB part

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
            dRq[IDX_R(egl,ell,m1,m2)] -= sqrt(a_sz) * ( sqrt(fnu[m1]*a_cib[ind2]) + sqrt(fnu[m2]*a_cib[ind1]) ) *
              corr_template[mell] * 2.0*M_PI/(dell*(dell+1.0));
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }
  return;
}

void sz_cib_cib_index_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err) {
  SZ_CIB_DEFS;

  SZ_CIB_INITS;

  parametric_zero_rq(egl,dRq);

  // NO SZ part

  // Add the CIB part
  if (do_cib == 1) {
    for (ell=egl->lmin;ell<=egl->lmax;ell++) {
      for(m1=0;m1<nfreq;m1++) {
        ind1 = ind_freq[m1];
        for (m2=m1;m2<nfreq;m2++){
          ind2 = ind_freq[m2];
          if ((ind1 >= 0) && (ind2 >= 0)) { // 100, 143 or 217 GHz
            dRq[IDX_R(egl,ell,m1,m2)] += r_cib[ind1*N_FREQ_CIB+ind2]*sqrt(a_cib[ind1]*a_cib[ind2])/d3000 * log((double)ell/(double)lnorm)*pow(((double)ell/(double)lnorm),cib_index);
            dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
          }
        }
      }
    }
  }

  // No CIB-SZ correlation part

  return;
}


// CIB a la GPE+Dunkley

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

void cib_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv,nrm;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;


  l_pivot = parametric_get_value(egl,"cib_l_pivot",err);
  forwardError(*err,__LINE__,);
  nrm = l_pivot*(l_pivot+1)/2/M_PI;

  index = parametric_get_value(egl,"cib_index",err);
  forwardError(*err,__LINE__,);

  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  stop = 0;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"cib_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = pow(ell/l_pivot,index);
          dRq[IDX_R(egl,ell,m1,m2)] = v/nrm;
          dRq[IDX_R(egl,ell,m2,m1)] = v/nrm;
        }
        break;
      }
    }
    if (stop==1) {
      return;
    }
  }      
  // error return
  parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

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

// PS a la GPE+Dunkley

void pointsource_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA,nrm;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  
  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = 1;
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

void pointsource_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrm;

  //l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  //forwardError(*err,__LINE__,);
  //
  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  stop = 0;
  v = 1/nrm;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          dRq[IDX_R(egl,ell,m1,m2)] = v;
          dRq[IDX_R(egl,ell,m2,m1)] = v;
        }
        break;
      }
    }
    if (stop==1) {
      return;
    }
  }      
  // error return
  parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

  return;
}

parametric *pointsource_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &pointsource_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ps_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  parametric_add_derivative_function(egl,"any",&pointsource_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}

// PS a la GPE+Dunkley

void pointsource_bydet_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA,nrm;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  
  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,(int)m2);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = 1;
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

void pointsource_bydet_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrm;

  //l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  //forwardError(*err,__LINE__,);
  //
  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  stop = 0;
  v = 1/nrm;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,m2);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          dRq[IDX_R(egl,ell,m1,m2)] = v;
          dRq[IDX_R(egl,ell,m2,m1)] = v;
        }
        break;
      }
    }
    if (stop==1) {
      return;
    }
  }      
  // error return
  parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

  return;
}

parametric *pointsource_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_bydet_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &pointsource_bydet_compute;
  forwardError(*err,__LINE__,NULL);
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ps_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)m1,(int)m2);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  parametric_add_derivative_function(egl,"any",&pointsource_bydet_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}


#define GPE_DUST_DEFS double gpe_dust_freqlist[3] = {100.,143.,217.}; \
  int nfreqs_gpe_dust = 3; \
  double A=28.2846, alpha=0.538389, B=657.4222, beta=2.96484, dellc=2029.09, gamma=1.68974, delta=42.1039; \
  double r_gpe_dust[3][3]; \
  double color = 3.16; \
  r_gpe_dust[0][0] = 0.; \
  r_gpe_dust[1][0] = 0.; \
  r_gpe_dust[0][1] = 0.; \
  r_gpe_dust[1][1] = 0.; \
  r_gpe_dust[2][0] = 0.; \
  r_gpe_dust[0][2] = 0.; \
  r_gpe_dust[2][1] = color; \
  r_gpe_dust[1][2] = color; \
  r_gpe_dust[2][2] = color*color;


// Dust component 'a la GPE'
void gpe_dust_compute(parametric *egl, double *Rq, error **err);

parametric *gpe_dust_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  GPE_DUST_DEFS;

  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  // Declare payload, allocate it and fill it

  parametric_check_freq(egl, gpe_dust_freqlist, nfreqs_gpe_dust,err);
  forwardError(*err,__LINE__,NULL);


  egl->payload = malloc_err(sizeof(template_payload),err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  // Get mapping between input frequencies and {143,217}
  payload->ind_freq = malloc_err(sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  for (m1=0;m1<egl->nfreq;m1++) {
    payload->ind_freq[m1]=-1; //Init 
    for (m2=0;m2<nfreqs_gpe_dust;m2++) {
      if (fabs(egl->freqlist[m1]-gpe_dust_freqlist[m2])<1e-6) {
        payload->ind_freq[m1]=m2;
      }
    }
  }
  

  // Allocate template and precompute it (value at 143 GHz)
  payload->template = malloc_err(sizeof(double)*(lmax-lmin+1),err);
  forwardError(*err,__LINE__,NULL);
  for (ell=lmin;ell<=lmax;ell++) {
    dell = (double) ell;
    payload->template[ell-lmin] = 1./(color*color-1.) * ( A*pow(100./dell,alpha) + B* pow(dell/1000.,beta)/pow(1.+pow(dell/dellc,gamma),delta) );
  //payload->template[ell-lmin] = 1./(color*color-1.) * ( A*exp(alpha*log(100./dell)) + B* exp(beta*log(dell/1000.))/exp(delta*log(1.+exp(gamma*log(dell/dellc)))) );
    payload->template[ell-lmin] *= 2.0*M_PI/dell/(dell+1.0); // Converts to C_ell
  }
  
  egl->eg_compute = &gpe_dust_compute;
  egl->eg_free = &parametric_template_payload_free;
  
  parametric_set_default(egl,"gpe_dust_norm",1.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gpe_dust_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gpe_dust_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  GPE_DUST_DEFS;

  nfreq = egl->nfreq;
  gpe_dust_norm = parametric_get_value(egl,"gpe_dust_norm",err);
  forwardError(*err,__LINE__,);
  payload = egl->payload;
  template = payload->template;

  ind_freq = payload->ind_freq;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<nfreq;m1++) {
      ind1 = ind_freq[m1];
      for (m2=m1;m2<nfreq;m2++) {
        ind2 = ind_freq[m2];
        if ((ind1 >=0) && (ind2 >=0)) { //100, 143 or 217
          Rq[IDX_R(egl,ell,m1,m2)] = gpe_dust_norm * r_gpe_dust[ind1][ind2] * template[ell-egl->lmin];
          Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
        } else {
          Rq[IDX_R(egl,ell,m1,m2)] = 0.0;
          Rq[IDX_R(egl,ell,m2,m1)] = 0.0;
        }
      }
    }
  }
  return;
}


CREATE_PARAMETRIC_FILE_INIT(radiogal,radiogal_init);
CREATE_PARAMETRIC_FILE_INIT(ir_clustered,ir_clustered_init);
CREATE_PARAMETRIC_FILE_INIT(ir_poisson,ir_poisson_init);
CREATE_PARAMETRIC_FILE_INIT(galametric,galactic_component_init);
CREATE_PARAMETRIC_FILE_INIT(poisson_tensor_bydet,poisson_tensor_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tensor_bydet,powerlaw_tensor_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(poisson_tensor,poisson_tensor_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tensor,powerlaw_tensor_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_triangle,powerlaw_triangle_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tanh,powerlaw_tanh_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_free_emissivity,powerlaw_free_emissivity_init);
CREATE_PARAMETRIC_FILE_INIT(cib,cib_init);
CREATE_PARAMETRIC_FILE_INIT(cibr,cibr_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz,sz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ksz,ksz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_cib,sz_cib_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_x,sz_x_init);
CREATE_PARAMETRIC_FILE_INIT(cib_x,cib_x_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz_cib_x,sz_cib_x_init);
CREATE_PARAMETRIC_FILE_INIT(pointsource,pointsource_init);
CREATE_PARAMETRIC_FILE_INIT(pointsource_bydet,pointsource_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(gpe_dust,gpe_dust_init);
