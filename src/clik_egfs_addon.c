#include "clik_helper.h"
#include "clik_egfs.h"

typedef struct {
  cmblkl* lkl;
  egfs *egfs_model;
  double *bars,*cls;
  double *rq;
  int lmin;
  } egfs_single;

double egfs_single_lkl(void* pelf,double *pars, error **err);
void egfs_single_free(void **pelf);
  
cmblkl* add_single_channel_egfs(cmblkl *base, int nvar, char **keyvars, int ndefaults, char** keys, char** values, double* cib_clustering,double *patchy_ksz, double *homogeneous_ksz,double *tsz, double *cib_decor_clust, double * cib_decor_poisson,error **err) {
  egfs_single *egfs_pay;
  cmblkl* res;
  int has_cl[20];
  int lmin,lmax,i,ii;
  char **xnames;
  
  lmin = base->ell[0];
  lmax = base->ell[base->nell-1];
  
  egfs_pay = malloc_err(sizeof(egfs_single),err);
  forwardError(*err,__LINE__,NULL);
  
  egfs_pay->egfs_model = egfs_init(nvar, keyvars, ndefaults, keys, values, lmin, lmax, cib_clustering,patchy_ksz, homogeneous_ksz,tsz, cib_decor_clust, cib_decor_poisson,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(egfs_pay->egfs_model->nfr!=1,-98765,"Bad number of channels fot egfs model (expected 1 got %d)",*err,__LINE__,NULL,egfs_pay->egfs_model->nfr);
  for(i=1;i<6;i++) {
    testErrorRet(base->offset_cl[i]!=-1,-345678,"Only pure T model at this time",*err,__LINE__,NULL);
  }
  
  egfs_pay->lkl = base;
  egfs_pay->lmin = lmin;
  egfs_pay->bars = malloc_err(sizeof(double)*(base->nbins+base->xdim),err);
  egfs_pay->cls = malloc_err(sizeof(double)*(base->nell),err);
  forwardError(*err,__LINE__,NULL);

  egfs_pay->rq = malloc_err(sizeof(double)*(lmax+1-lmin),err);
  forwardError(*err,__LINE__,NULL);
  
  memset(has_cl,0,sizeof(int)*20);
  has_cl[0]=1;
  
  res = init_cmblkl(egfs_pay, egfs_single_lkl, 
                      egfs_single_free,
                      base->nell,base->ell,has_cl,-1,base->unit,base->wl,0,
                      base->bins,base->nbins, base->xdim+nvar,err);
  forwardError(*err,__LINE__,NULL);

  testErrorRet(base->xdim!=0 && base->xnames==NULL,-5864,"XNames unset",*err,__LINE__,NULL);
  
  xnames = malloc_err(sizeof(char*)*(base->xdim+nvar),err);
  forwardError(*err,__LINE__,NULL);
  
  for(ii=0;ii<base->xdim;ii++) {
    xnames[ii] = base->xnames[ii];
  }
  for(ii=0;ii<nvar;ii++) {
    xnames[ii+base->xdim] = keyvars[ii];
  }
  
  cmblkl_set_names(res, xnames,err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames);
  
  return res;
}

double egfs_single_lkl(void* pelf,double *pars, error **err) {
  egfs_single *egfs_pay;
  double *wl,*wl0,one;
  int inc,il;
  double res;
  
  egfs_pay = pelf;
  egfs_compute(egfs_pay->egfs_model, &(pars[egfs_pay->lkl->nbins+egfs_pay->lkl->xdim]), egfs_pay->rq, NULL, err);
  forwardError(*err,__LINE__,0);
  
  // apply wl and binning
  one=1;
  if (egfs_pay->lkl->wl==NULL) {
    wl0 = &one;
    inc = 0;
  } else {
    wl0 = egfs_pay->lkl->wl;
    inc = 1;
  }
  
  wl=wl0;
  for(il=0;il<egfs_pay->lkl->nell;il++) {
    egfs_pay->cls[il] = egfs_pay->rq[egfs_pay->lkl->ell[il]-egfs_pay->lmin] * *wl * egfs_pay->lkl->unit;
    wl+=inc;
  }
  
  memcpy(egfs_pay->bars,pars,sizeof(double)*(egfs_pay->lkl->nbins+egfs_pay->lkl->xdim));
  
  // apply binning if needed
  if (egfs_pay->lkl->bins!=NULL) {
    char trans;
    int npar;
    double done,dzero;
    int one;
    int ndim;
    
    trans='T';
    npar = egfs_pay->lkl->nbins;
    one = 1;
    done = 1.;
    dzero = 0;
    ndim = egfs_pay->lkl->ndim;
    
    
    //_DEBUGHERE_("cls[0]=%g pls[0]=%g bins[0]=%g",cls[0],llkl->pls[0],llkl->bins[0]);
    
    dgemv(&trans, &ndim, &npar, &done, egfs_pay->lkl->bins, &ndim, egfs_pay->cls, &one, &done, egfs_pay->bars, &one);
    //_DEBUGHERE_("cls[0]=%g pls[0]=%g ",cls[0],llkl->pls[0]);
  } else {
    for(il=0;il<egfs_pay->lkl->nell;il++) {
      egfs_pay->bars[il] += egfs_pay->cls[il];
    }
  }
  
  res = egfs_pay->lkl->lkl_func(egfs_pay->lkl->lkl_data,egfs_pay->bars,err);
  forwardError(*err,__LINE__,0);
  
  return res;
}

void egfs_single_free(void **pelf) {
  egfs_single *egfs_pay;
  
  egfs_pay = *pelf;
  free(egfs_pay->bars);
  free(egfs_pay->cls);
  free(egfs_pay->rq);
  
  egfs_free(&(egfs_pay->egfs_model));
  free_cmblkl(&(egfs_pay->lkl));
  
  free(egfs_pay);
  
  *pelf=NULL;
}
