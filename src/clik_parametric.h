#include "errorlist.h"
#include "io.h"
#include <math.h>

#ifndef CLIK_PARAMETRIC
#define CLIK_PARAMETRIC


#define pfcharsize 256

typedef char pfchar[pfcharsize];

typedef struct {
  pfchar *key;
  pfchar *value;
  int nkey;
  int nmax;
} pflist;

pflist* pflist_init(error **err);
void pflist_add_item(pflist* pf,int nit, char** key, char **value,error **err);
void pflist_free(void **ppf);
char* pflist_get_value(pflist* pf, char* key,char* safeguard, error **err);
int pflist_key_index(pflist *pf, char *key, error **err);
void pflist_remove_item(pflist* pf, int index,error **err);
long pflist_get_int_value(pflist *pf, char *key,long* safeguard, error **err);
double pflist_get_double_value(pflist *pf, char *key,double *safeguard, error **err);

typedef void (exg_compute)(void* exg, double *Rq, double*dRq, error **err);
typedef void (exg_freepayload)(void** data);

typedef struct {
  void *payload;
  pflist *pf;
  int *freqlist;
  int *det2freq;
  int ndet;
  int nfreq;
  exg_compute *eg_compute;
  exg_freepayload *eg_free;
  double *sRq, *sdRq;
  int nvar;
  int ndef;
  int lmin,lmax;
  pfchar *varkey;
  pflist *default_settings;
  int dnofail;
} parametric;

parametric *parametric_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void parametric_free(void** pegl);
void parametric_compute(parametric *egl, double *pars, double* Rq, double *dRq, error **err);
double parametric_get_default(parametric* egl,char *key, error **err);
double parametric_get_value(parametric *egl, char *key, error **err);
void parametric_dnofail(parametric* egl, int vl);


parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void powerlaw_compute(void* exg, double *Rq, double*dRq, error **err);
parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void powerlaw_free_emissivity_free(void **pp);
void powerlaw_free_emissivity_compute(void* exg, double *Rq, double*dRq, error **err);

parametric *radiogal_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void radiogal_compute(void* exg, double *Rq, double* dRq, error **err);
double dBdT(double nu, double nu0);
void radiogal_free(void **pp);
double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0);
double d_dust_spectrum_d_beta_dust(double nu, double T_dust, double beta_dust, double nu0);
double d_dust_spectrum_d_T_dust(double nu, double T_dust, double beta_dust, double nu0);
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0);
double d_non_thermal_spectrum_d_alpha_non_thermal(double nu, double alpha_non_thermal, double nu0);
parametric *galactic_component_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void galactic_component_compute(void* exg, double *Rq, double *dRq, error **err);
void galactic_component_free(void **pp);
parametric *ir_poisson_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void ir_poisson_compute(void* exg, double *Rq, double *dRq, error **err);
void ir_poisson_free(void **pp);
parametric *ir_clustered_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void ir_clustered_compute(void* exg, double *Rq, double *dRq, error **err);
void ir_clustered_free(void **pp);

#endif
