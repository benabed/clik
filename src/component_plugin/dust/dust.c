// group the dust stuff here
#include "dust.h"



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



// Non-thermal emissivity normalized to 1 at nu0 (e.g. synchrotron, free-free, etc.)
// The spectral index is called alpha here, to distinguish from the gray body case
// Expressed in dT (CMB)
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0) {

  //return (pow(nu/nu0,alpha_non_thermal)/dBdT(nu,nu0));
  return (exp(alpha_non_thermal*log(nu/nu0))/dBdT(nu,nu0));
}

// TT

void hgal_compute(parametric *egl, double *Rq, error **err) {
  double *A;
  int nfreq,m1,m2,ell;
  double a,b,v,l_pivot,lA;
  pfchar name;
  int mell;

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"hgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"hgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  b = parametric_get_value(egl,"hgal_beta",err);
  forwardError(*err,__LINE__,);
  l_pivot = parametric_get_value(egl,"hgal_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = exp(b*(ell-l_pivot))/l_pivot/(l_pivot+1)*2*M_PI;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }
}

parametric *hgal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &hgal_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"hgal_l_pivot",250,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"hgal_beta",-0.00849,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"hgal_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"hgal_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,50,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  
  return egl;
}
CREATE_PARAMETRIC_FILE_INIT(hgal,hgal_init);

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
  forwardError(*err,__LINE__,NULL);

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
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,NULL,type);
    // return ?
  }

  return egl;
}

// P definitions
typedef struct {
  int n;
  int *m1,*m2;
  double *nrm;
  char **nrm_names;
} pw_XX_payload;

#define TE_KIND 1
#define EE_KIND 2
pw_XX_payload*  init_pw_XX_payload(int kind,int nT,int nP,error **err);

void pw_XX_free(void **PP);

//TE
void gal_TE_compute(parametric* exg, double *Rq, error **err);

parametric *gal_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &gal_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,err);
  forwardError(*err,__LINE__,NULL);


  parametric_set_default(egl,"gal_TE_norm",1,err);
  forwardError(*err,__LINE__,NULL);
    
  parametric_set_default(egl,"gal_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_TE_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_beta_dust",1.8,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_TE_T_dust",18,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gal_TE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v,norm;
  double prm_nu0,beta_dust,T_dust;

  payload = egl->payload;
  mx = egl->nfreq_T + egl->nfreq_P;
  

  prm_nu0 = parametric_get_value(egl,"gal_TE_nu0",err);
  forwardError(*err,__LINE__,);

  beta_dust = parametric_get_value(egl,"gal_TE_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_TE_T_dust",err);
  forwardError(*err,__LINE__,);

  for(i=0;i<mx;i++) {
    payload->nrm[i] = dust_spectrum((double)egl->freqlist[i],T_dust,beta_dust,prm_nu0);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"gal_TE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"gal_TE_index",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_TE_norm",err);
  forwardError(*err,__LINE__,);
  

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index) * norm;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1] * A[m2] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

//EE
void gal_EE_compute(parametric* exg, double *Rq, error **err);

parametric *gal_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta;
  double prm_nu0,beta_dust,T_dust;

  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &gal_EE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(EE_KIND,egl->nfreq_T,egl->nfreq_P,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_norm",1,err);
  forwardError(*err,__LINE__,NULL);
    
  parametric_set_default(egl,"gal_EE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"gal_EE_nu0",PRM_NU0,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_beta_dust",1.8,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_EE_T_dust",18,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gal_EE_compute(parametric* egl, double *Rq, error **err) {
  pw_XX_payload *payload;
  int i,mx;
  int ell,m1,m2;
  double *A;
  double l_pivot, index,v,norm;
  double prm_nu0,beta_dust,T_dust;

  payload = egl->payload;
  mx = egl->nfreq_P;
  

  prm_nu0 = parametric_get_value(egl,"gal_EE_nu0",err);
  forwardError(*err,__LINE__,);

  beta_dust = parametric_get_value(egl,"gal_EE_beta_dust",err);
  forwardError(*err,__LINE__,);

  T_dust = parametric_get_value(egl,"gal_EE_T_dust",err);
  forwardError(*err,__LINE__,);

  for(i=0;i<mx;i++) {
    payload->nrm[i] = dust_spectrum((double)egl->freqlist[i+egl->nfreq_T],T_dust,beta_dust,prm_nu0);
    forwardError(*err,__LINE__,);
  }
  A = payload->nrm;

  l_pivot = parametric_get_value(egl,"gal_EE_l_pivot",err);
  forwardError(*err,__LINE__,);

  index = parametric_get_value(egl,"gal_EE_index",err);
  forwardError(*err,__LINE__,);

  norm = parametric_get_value(egl,"gal_EE_norm",err);
  forwardError(*err,__LINE__,);
  

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index) * norm;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      Rq[IDX_R(egl,ell,m1,m2)] = A[m1-egl->nfreq_T] * A[m2-egl->nfreq_T] * v;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
    }
  }
}

// definitions
CREATE_PARAMETRIC_FILE_INIT(galametric,galactic_component_init);
CREATE_PARAMETRIC_POL_FILE_INIT(gal_TE,gal_TE_init);
CREATE_PARAMETRIC_POL_FILE_INIT(gal_EE,gal_EE_init);

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

CREATE_PARAMETRIC_FILE_INIT(gpe_dust,gpe_dust_init);

