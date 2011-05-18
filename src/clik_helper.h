#include "clik.h"
#include "hdf5.h"
#include "lklbs.h"
#include "aplowly.h"
#include "fowly.h"
#include "smica.h"
#include <dlfcn.h>

#ifndef _CLIK_HLP_
#define _CLIK_HLP_

#define _dealwitherr error *lerr,**err; if(_err==NULL) {lerr=NULL;err=&lerr;} else {err=_err;}

#define _forwardError(A,B,C) if(_err!=NULL) {forwardError(A,B,C);} else {quitOnError(A,B,stderr);}
#define _testErrorRetVA(A,B,C,D,E,F,...) if(_err!=NULL) {testErrorRetVA(A,B,C,D,E,F,__VA_ARGS__);} else {testErrorExitVA(A,B,C,D,E,__VA_ARGS__);}

// hdf helpers
double* hdf5_double_datarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err);
int* hdf5_int_datarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err);
double* hdf5_double_attarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err);

// get environ parameter
int clik_getenviron_integer(char* name, int sfg, error **err);
double clik_getenviron_real(char* name, double sfg, error **err);
char* clik_getenviron_string(char* name, char* sfg, error **err);
int clik_getenviron_numthread(char* name, int sfg, error **err);

// init lkls
cmblkl * clik_lklobject_init(hid_t group_id,char* cur_lkl,error **err);
typedef cmblkl* clik_lkl_init_func(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err);
#endif