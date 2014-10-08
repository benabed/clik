#include "clik_helper.h"

typedef struct {
  double *datavector;
  double *bin;
  double *siginv;
  double *temp,*temp2;
  int n_data;
  int n_input;
} dataset_struct;


double dataset_lkl(void *vds, double* pars,error **err) {
  dataset_struct *ds;
  int i;
  int j;
  double lkl;

  vds = ds;
  for(i=0;i<ds->n_data;i++) {
    ds->temp[i] = -ds->datavector[i];
    for(j=0;j<ds->n_input;i++) {
      ds->temp[i] += ds->bin[i*ds->n_input+j] * pars[j];
    }
  }

  lkl = 0;
  for(i=0;i<ds->n_data;i++) {
    ds->temp2[i] = 0;
    for(j=0;i<ds->n_data;j++) {
      ds->temp2[i] += ds->siginv[i*ds->n_data+j] * ds->temp[j];
    }
    lkl += ds->temp2[i]*ds->temp[i];
  }

  return -.5*lkl;
}

void* dataset_init(int n_data, double *data, int n_input, double* bin, double* siginv, error **err ) {
  dataset_struct *ds;

  ds = malloc_err(sizeof(dataset_struct),err);
  forwardError(*err,__LINE__,NULL);

  ds->n_data = n_data;
  ds->n_input = n_input;

  ds->bin = malloc_err(sizeof(double)*n_data*n_input,err);
  forwardError(*err,__LINE__,NULL);

  ds->siginv = malloc_err(sizeof(double)*n_data*n_data,err);
  forwardError(*err,__LINE__,NULL);

  ds->datavector = malloc_err(sizeof(double)*n_data,err);
  forwardError(*err,__LINE__,NULL);

  ds->temp = malloc_err(sizeof(double)*n_data,err);
  forwardError(*err,__LINE__,NULL);

  ds->temp2 = malloc_err(sizeof(double)*n_data,err);
  forwardError(*err,__LINE__,NULL);

  return ds;
}

void dataset_free(void **vvds) {
  dataset_struct *ds;

  ds = *vvds;

  free(ds->datavector);
  free(ds->bin);
  free(ds->siginv);
  free(ds->temp);
  free(ds->temp2);
  free(ds);
  *vvds = NULL;
}

