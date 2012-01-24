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
} parametric;

parametric *parametric_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void parametric_free(void** pegl);
void parametric_compute(parametric *egl, double *pars, double* Rq, double *dRq, error **err);

parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void powerlaw_compute(void* exg, double *Rq, double*dRq, error **err);
parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err);
void powerlaw_free_emissivity_free(void **pp);
void powerlaw_free_emissivity_compute(void* exg, double *Rq, double*dRq, error **err);

#endif