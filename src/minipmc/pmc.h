/*
 *  minipmc.h
 *  simplified pmc support from pmclib
 */


#ifndef __MC_DIST
#define __MC_DIST

#include "errorlist.h"
#include "io.h"

#include <dlfcn.h>

struct _pmc_simu_struct_;
struct _distribution_struct_;

typedef double posterior_log_pdf_func(void *, double *, error **);
typedef posterior_log_pdf_func log_pdf_func;

typedef void retrieve_ded_func(const void *, double *, error **);
typedef void (posterior_log_free)(void**);
typedef posterior_log_free free_func;

typedef long simulate_func(struct _pmc_simu_struct_ *, void *, void *, void *, error **);
typedef void filter_func(struct _pmc_simu_struct_* , void *, error **);
typedef void update_func(void *, struct _pmc_simu_struct_ *, error **);
typedef void* mpi_exchange_func(void *, error **);
typedef double first_derivative_func(void*, int , const double*, error **err);
typedef double second_derivative_func(void*, int ,int, const double*, error **err);
 
typedef char _char_name[1024];

typedef struct _distribution_struct_ {
  int ndim, n_ded,ndef;
  double *pars;
  int * def;
  void* data;
  posterior_log_pdf_func* log_pdf;
  retrieve_ded_func* retrieve;
  posterior_log_free* free;
  simulate_func *simulate;
  mpi_exchange_func *broadcast_mpi;
  first_derivative_func *f_der;
  second_derivative_func *d_der;
  
  void* dlhandle;
  _char_name *name;
   
} distribution;

distribution* init_distribution_full(int ndim,
                                     void* data, 
                                     posterior_log_pdf_func* log_pdf,
                                     posterior_log_free* freef,
                                     simulate_func *simulate,
                                     int nded,
                                     retrieve_ded_func* retrieve,
                                     error **err);

distribution* init_simple_distribution(int ndim,
                               void* data, 
                               posterior_log_pdf_func* log_pdf,
                               posterior_log_free* freef,
                               error **err);

distribution* init_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                simulate_func *simulate,
                                error **err);
void free_distribution(distribution **pdist) ;

int distribution_get_name(distribution *dist,char* name,error **err);
int* distribution_get_names(distribution *dist,int nname,char** name,int includeded,error **err);
void distribution_set_names(distribution *dist,char** name, error **err);

void distribution_set_broadcast(distribution* dist, mpi_exchange_func* broadcast, error **err);


double distribution_lkl(void* pdist, const double* pars, error **err);
#define distribution_log_pdf distribution_lkl

void distribution_retrieve(const void* pdist, double* pars, error **err);

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err);
void distribution_set_default_name(distribution *dist, int ndef, char** idef, double* vdef,error **err);

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err);


#define pmc_base       -6000
#define pmc_allocate   -1 + pmc_base
#define pmc_serialize  -2 + pmc_base
#define pmc_outOfBound -3 + pmc_base
#define pmc_badComm    -4 + pmc_base
#define pmc_negWeight  -5 + pmc_base
#define pmc_cholesky   -6 + pmc_base
#define pmc_negative   -7 + pmc_base
#define pmc_undef      -8 + pmc_base
#define pmc_file       -9 + pmc_base
#define pmc_io        -10 + pmc_base
#define pmc_tooManySteps -11 + pmc_base
#define pmc_dimension -12 + pmc_base
#define pmc_type      -13 + pmc_base
#define pmc_negHatCl  -14 + pmc_base
#define pmc_infnan    -15 + pmc_base
#define pmc_incompat  -16 + pmc_base
#define pmc_nosamplep -17 + pmc_base
#define pmc_sort      -18 + pmc_base
#define pmc_infinite  -19 + pmc_base
#define pmc_isLog     -20 + pmc_base

#define dist_base       -6700
#define dist_undef       -1 + dist_base 
#define dist_type        -2 + dist_base 

#define mv_base       -700
#define mv_allocate   -1 + mv_base
#define mv_serialize  -2 + mv_base
#define mv_outOfBound -3 + mv_base
#define mv_badComm    -4 + mv_base
#define mv_negWeight  -5 + mv_base
#define mv_cholesky   -6 + mv_base
#define mv_negative   -7 + mv_base
#define mv_undef      -8 + mv_base
#define mv_file       -9 + mv_base
#define mv_io        -10 + mv_base
#define mv_tooManySteps -11 + mv_base
#define mv_dimension -12 + mv_base
#define mv_type      -13 + mv_base
#define mv_negHatCl  -14 + mv_base

#define PI         3.141592653589793
#define LOGSQRT2PI 0.918938533204673

#define MALLOC_IF_NEEDED(ret,dest,size,err) { \
  if (dest==NULL) { \
    ret= (double*) malloc_err(size,err); \
  } else { \
    ret=dest; \
  } \
}

#endif
