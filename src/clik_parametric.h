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
  double *dvalue;
  int nkey;
  int nmax;
  int ncommon;
} pflist;

pflist* pflist_init(error **err);
void pflist_add_item(pflist* pf,int nit, char** key, char **value,error **err);
void pflist_free(void **ppf);
char* pflist_get_value(pflist* pf, char* key,char* safeguard, error **err);
int pflist_key_index(pflist *pf, char *key, error **err);
void pflist_remove_item(pflist* pf, int index,error **err);
long pflist_get_int_value(pflist *pf, char *key,long* safeguard, error **err);
double pflist_get_double_value(pflist *pf, char *key,double *safeguard, error **err);
void pflist_compute_ncommon(pflist *pf, error **err);

struct parametric_struct;

typedef void (exg_compute)(struct parametric_struct* exg, double *Rq,error **err);
typedef void (exg_deriv)(struct parametric_struct* exg, int iv,double *Rq, double*dRq, error **err);
typedef void (exg_compute_and_deriv)(struct parametric_struct* exg, double *Rq, double*dRq, error **err);
typedef void (exg_freepayload)(void** data);

typedef struct parametric_struct {
  void *payload;
  pflist *pf;
  double *freqlist;
  double *detlist;
  int *det2freq;
  int ndet;
  int nfreq;
  exg_compute *eg_compute;
  exg_compute_and_deriv *eg_compute_and_deriv;
  exg_deriv **eg_deriv;
  exg_deriv *eg_deriv_any;
  pfchar *deriv_key;
  int nderiv;
  exg_freepayload *eg_free;
  double *sRq, *sdRq;
  int nvar;
  int ndef;
  int lmin,lmax;
  pfchar *varkey;
  pflist *default_settings;
  int dnofail;
  double *color;
  double l_pivot;
  pfchar tensor_norm_template;
  int tensor_norm_template_len;
} parametric;

parametric *parametric_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
parametric *parametric_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void parametric_free(void** pegl);
void parametric_compute(parametric *egl, double *pars, double* Rq, double *dRq, error **err);
void parametric_compute_loop(parametric * egl, double* Rq, double *dRq, error **err);
double parametric_get_default(parametric* egl,char *key, error **err);
double parametric_get_value(parametric *egl, char *key, error **err);
void parametric_dnofail(parametric* egl, int vl);
void parametric_set_default(parametric* egl,char *key, double value,error **err);
void parametric_add_derivative_function(parametric *egl, char* vk, exg_deriv* fnc,error **err);
void parametric_end_derivative_loop(parametric *egl,double* dRq, char* varkey, error **err);
void parametric_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);
void parametric_index_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);

void parametric_set_color(parametric *egl,double *color, error **err);

void parametric_check_freq(parametric *egl, double* frq, int nfreq, error **err);

double dBdT(double nu, double nu0);

parametric *powerlaw_init(int ndet, double* detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
parametric *powerlaw_free_emissivity_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);

#define PTR_DER(egl,iv,dRq) (&(dRq[iv*(egl->lmax+1-egl->lmin)*egl->nfreq*egl->nfreq]))
#define IDX_R(egl,ell,m1,m2) ((ell-egl->lmin)*egl->nfreq*egl->nfreq + m1*egl->nfreq + m2)

typedef struct {
  double *template;
  int *ind_freq;
} template_payload;
void parametric_template_payload_init(parametric *egl, double *template, int template_size, double *freqlist, int nfreq, error **err);
void parametric_template_payload_free(void **pp);

void parametric_simple_payload_free(void **pp);

void parametric_tensor_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);
void parametric_tensor_fill(parametric *egl,double *A,error **err);

#endif
