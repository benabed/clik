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


typedef struct {
  int nbins,lmax,nlt;
  double *bins;
  double *cors, *cor0;
  double *siginv;
  double *buf,*p_hat;
  int renorm,has_calib;
  double *cl_fid;
} dts_lensing;

double dts_lensing_lkl(void *vds,double *pars, error **err) {
  dts_lensing *dt;
  int il,ib,it;
  double b_ib;
  int cut;
  double lkl,bv;
  double calib;
  double *fpars;

  dt = vds;

  calib = 1;
  if (dt->has_calib==1) {
    calib = pars[dt->nlt];
    calib = 1./(calib*calib);  
  }
  
  fpars = pars;
  if (dt->renorm == 0) {
    fpars = dt->cl_fid;
  }

  for(ib=0;ib<dt->nbins;ib++) {
    b_ib = 0;
    
    cut=0;
    for(il=0;il<dt->lmax+1;il++) {
      bv = dt->bins[ib*(dt->lmax+1) + il];
      if (bv==0) {
        if (cut==0) {
          continue;
        } else {
          break;
        }
      }
      b_ib += (bv * pars[il]) * calib;
      cut=1;
    }

    if (dt->cor0!=NULL) {
      b_ib -= dt->cor0[ib];
    }
    
    for(il=0; il<dt->nlt;il++) {
      b_ib += (dt->cors[ib*dt->nlt + il] * fpars[il])*calib;
    }

    dt->buf[ib] = dt->p_hat[ib] - b_ib;
  }
  
  lkl = 0;

  for(ib=0;ib<dt->nbins;ib++) {
    lkl += dt->siginv[ib*dt->nbins+ib] * dt->buf[ib];
    for(it=ib+1;it<dt->nbins;it++) {
      lkl += 2*dt->siginv[ib*dt->nbins+it] * dt->buf[it];
    }
  }

  return -.5*lkl;

}

void* dts_lensing_init(int lmax,int *hascl, int nbins, double* bins, double *p_hat, double *cors, double *cor0, double* siginv,double *cl_fid, int renorm ,int has_calib, error **err ) {
  dts_lensing *ds;
  int cli;

  ds = malloc_err(sizeof(dataset_struct),err);
  forwardError(*err,__LINE__,NULL);

  ds->nbins = nbins;
  ds->lmax = lmax;
  ds->renorm = renorm;
  ds->has_calib = has_calib;
  
  ds->bins = malloc_err(sizeof(double)*nbins*(lmax+1),err);
  forwardError(*err,__LINE__,NULL);

  memcpy(ds->bins,bins,sizeof(double)*nbins*(lmax+1));

  ds->p_hat = malloc_err(sizeof(double)*nbins,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(ds->p_hat,p_hat,sizeof(double)*nbins);

  ds->siginv = malloc_err(sizeof(double)*nbins*nbins,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(ds->siginv,siginv,sizeof(double)*nbins*nbins);

  ds->buf = malloc_err(sizeof(double)*nbins,err);
  forwardError(*err,__LINE__,NULL);


  ds->nlt = 0;
  ds->cors = NULL;
  tlt = lmax + 1;

  if (cors!=NULL) {
    ds->nlt = lmax+1;
    for(cli=0;cli<6;cli++) {
      if (hascl[cli]!=0) {
        ds->nlt += lmax+1;
      }
    }
    ds->cors = malloc_err(sizeof(double)*nbins*ds->nlt,err);
    forwardError(*err,__LINE__,NULL);

    memcpy(ds->cors,cors,sizeof(double)*nbins*ds->nlt);

    tlt = ds->nlt;
  } 

  ds->cl_fid = malloc_err(sizeof(double)*ds->tlt,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(ds->cl_fid,cl_fid,sizeof(double)*ds->tlt);


  ds->cor0 = NULL;
  if (cor0!=NULL) {
    ds->cor0 = malloc_err(sizeof(double)*nbins,err);
    forwardError(*err,__LINE__,NULL);

    memcpy(ds->cor0,cor0,sizeof(double)*nbins);
  }
  
  return ds;
}


void dts_lensing_free(void **vvds) {
  dts_lensing *ds;

  ds = *vvds;

  free(ds->p_hat);
  free(ds->bins);
  free(ds->siginv);
  free(ds->buf);
  if (ds->cor0!=NULL) {
    free(ds->cor0);
  }
  if (ds->cors!=NULL) {
    free(ds->cors);
  }

  free(ds);
  *vvds = NULL;
}
