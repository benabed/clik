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
  forwardError(*err,__LINE__,);

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

  return (pow(nu/nu0,alpha_non_thermal)/dBdT(nu,nu0));
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
      a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,PRM_NU0);
    } else {
      a[m1] = non_thermal_spectrum((double)egl->freqlist[m1],alpha_non_thermal,PRM_NU0);
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

  // dR/dbeta_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,PRM_NU0);
    b[m1] = d_dust_spectrum_d_beta_dust((double)egl->freqlist[m1],T_dust,beta_dust,PRM_NU0);
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = a[m1]*b[m2]+a[m2]*b[m1];
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

  // dR/dT_dust
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,PRM_NU0);
    b[m1] = d_dust_spectrum_d_T_dust((double)egl->freqlist[m1],T_dust,beta_dust,PRM_NU0);
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = a[m1]*b[m2]+a[m2]*b[m1];
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

  // dR/dalpha_non_thermal
  // Get vector emissivity derivative
  for (m1=0;m1<nfreq;m1++) {
    a[m1] = non_thermal_spectrum((double)egl->freqlist[m1],alpha_non_thermal,PRM_NU0);
    b[m1] = d_non_thermal_spectrum_d_alpha_non_thermal((double)egl->freqlist[m1],alpha_non_thermal,PRM_NU0);
    for (m2=m1;m2<nfreq;m2++) {
      B[m1*nfreq+m2] = a[m1]*b[m2]+a[m2]*b[m1];
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
  egl->eg_compute = &galactic_component_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*2*egl->nfreq*(egl->nfreq+1),err);
  forwardError(*err,__LINE__,NULL);
  
  // uK^2 at l=500, nu=143 GHz;
  parametric_set_default(egl,"gal_norm",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_norm",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_norm",parametric_norm_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_index",0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_declare_mandatory(egl,"gal_index",err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"gal_index",parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,);

  if (strcmp(pt,"dust")==0) {
    parametric_set_default(egl,"gal_beta_dust",1.8,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_beta_dust",gal_beta_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

    parametric_set_default(egl,"gal_T_dust",18,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_T_dust",gal_T_dust_derivative,err);  
    forwardError(*err,__LINE__,NULL);

  } else if (strcmp(pt,"non_thermal")==0) {
    //Intensity, = -3.0 in RJ
    parametric_set_default(egl,"gal_alpha_non_thermal",-1,err);
    forwardError(*err,__LINE__,NULL);
    parametric_add_derivative_function(egl,"gal_alpha_non_thermal",gal_alpha_non_thermal_derivative,err);  
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

CREATE_PARAMETRIC_FILE_INIT(radiogal,radiogal_init);
CREATE_PARAMETRIC_FILE_INIT(ir_clustered,ir_clustered_init);
CREATE_PARAMETRIC_FILE_INIT(ir_poisson,ir_poisson_init);
CREATE_PARAMETRIC_FILE_INIT(galametric,galactic_component_init);
CREATE_PARAMETRIC_FILE_INIT(poisson_tensor_bydet,poisson_tensor_bydet_init);
CREATE_PARAMETRIC_FILE_INIT(powerlaw_tensor_bydet,powerlaw_tensor_bydet_init);
