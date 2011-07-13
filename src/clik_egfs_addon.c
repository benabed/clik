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

egfs* clik_egfs_init_hdf5(hid_t group_id,char* cur_lkl,error **err) {
  egfs *egfs_model;
  int nvar,ndef;
  herr_t hstat;
  int dz,i,lmin,lmax;
  char *keyvartable,*deftable,*valtable;
  char **keyvars,**keys,**values;
  double* cib_clustering,*patchy_ksz, *homogeneous_ksz, *tsz, *cib_decor_clust,*cib_decor_poisson;
  
  // get variable names
  hstat = H5LTget_attribute_int( group_id, ".", "ndim",  &nvar);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read ndim in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  keyvartable = hdf5_char_attarray(group_id,cur_lkl,"keys",&dz, err);
  forwardError(*err,__LINE__,NULL);
  
  keyvars = malloc_err(sizeof(char*)*nvar,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<nvar;i++) {
    keyvars[i] = &(keyvartable[i*256]);
  }  
  
  // get defaults
  hstat = H5LTget_attribute_int( group_id, ".", "ndef",  &ndef);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read ndef in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  deftable = hdf5_char_attarray(group_id,cur_lkl,"defaults",&dz, err);
  forwardError(*err,__LINE__,NULL);
  valtable = hdf5_char_attarray(group_id,cur_lkl,"defaults",&dz, err);
  forwardError(*err,__LINE__,NULL);
  
  keys = malloc_err(sizeof(char*)*ndef,err);
  forwardError(*err,__LINE__,NULL);
  values = malloc_err(sizeof(char*)*ndef,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<ndef;i++) {
    keys[i] = &(deftable[i*256]);
    values[i] = &(valtable[i*256]);
  }  
  
  // get lmin and lmax
  hstat = H5LTget_attribute_int( group_id, ".", "lmin",  &lmin);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax",  &lmax);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  // get data
  
  cib_clustering = hdf5_double_datarray(group_id,cur_lkl,"cib_clustering_template",&dz, err);
  forwardError(*err,__LINE__,NULL);

  patchy_ksz = hdf5_double_datarray(group_id,cur_lkl,"patchy_ksz_template",&dz, err);
  forwardError(*err,__LINE__,NULL);
  
  homogeneous_ksz = hdf5_double_datarray(group_id,cur_lkl,"homogeneous_ksz_template",&dz, err);
  forwardError(*err,__LINE__,NULL);
  
  tsz = hdf5_double_datarray(group_id,cur_lkl,"tsz_template",&dz, err);
  forwardError(*err,__LINE__,NULL);
  
  cib_decor_clust = NULL;
  hstat = H5LTfind_attribute (group_id, "cib_decor_clustering");
  if (hstat==1) {
    cib_decor_clust = hdf5_double_datarray(group_id,cur_lkl,"cib_decor_clustering",&dz, err);
    forwardError(*err,__LINE__,NULL);
  }
  
  cib_decor_poisson = NULL;
  hstat = H5LTfind_attribute (group_id, "cib_decor_poisson");
  if (hstat==1) {
    cib_decor_poisson = hdf5_double_datarray(group_id,cur_lkl,"cib_decor_poisson",&dz, err);
    forwardError(*err,__LINE__,NULL);
  }

  egfs_model = egfs_init( nvar, keyvars, ndef, keys, values, 
                         lmin, lmax, cib_clustering,patchy_ksz, 
                         homogeneous_ksz,tsz,
                         cib_decor_clust, cib_decor_poisson,err);
  forwardError(*err,__LINE__,NULL);

  free(keyvartable);
  free(keyvars);
  free(deftable);
  free(keys);
  free(valtable);
  free(values);
  free(cib_clustering);
  free(patchy_ksz);
  free(homogeneous_ksz);
  free(tsz);
  if (cib_decor_clust!=NULL) {
    free(cib_decor_clust);
  }        
  if (cib_decor_poisson!=NULL) {
    free(cib_decor_poisson);
  }        
  
  return egfs_model;     
}

cmblkl *clik_addon_egfs_single_init(cmblkl *base, hid_t group_id,char* cur_lkl,error **err) {
  egfs *egfs_model;
  egfs_single *egfs_pay;
  cmblkl* res;
  int has_cl[20];
  int lmin,lmax,i,ii;
  char **xnames;
  
  egfs_model = clik_egfs_init_hdf5(group_id,cur_lkl,err);
  forwardError(*err,__LINE__,NULL);
  
  lmin = base->ell[0];
  lmax = base->ell[base->nell-1];
  testErrorRetVA(egfs_model->nfr!=1,-98765,"Bad number of channels fot egfs model (expected 1 got %d)",*err,__LINE__,NULL,egfs_model->nfr);
  testErrorRetVA(egfs_model->lmin!=lmin,-98765,"Bad lmin (expected %d got %d)",*err,__LINE__,NULL,lmin,egfs_model->lmin);
  testErrorRetVA(egfs_model->lmax!=lmax,-98765,"Bad lmin (expected %d got %d)",*err,__LINE__,NULL,lmax,egfs_model->lmax);
  for(i=1;i<6;i++) {
    testErrorRet(base->offset_cl[i]!=-1,-345678,"Only pure T model at this time",*err,__LINE__,NULL);
  }
  
  egfs_pay = malloc_err(sizeof(egfs_single),err);
  forwardError(*err,__LINE__,NULL);
  
  egfs_pay->egfs_model = egfs_model;  
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
                      base->bins,base->nbins, base->xdim+egfs_model->np,err);
  forwardError(*err,__LINE__,NULL);

  testErrorRet(base->xdim!=0 && base->xnames==NULL,-5864,"XNames unset",*err,__LINE__,NULL);
  
  xnames = malloc_err(sizeof(char*)*(base->xdim+egfs_model->np),err);
  forwardError(*err,__LINE__,NULL);
  
  for(ii=0;ii<base->xdim;ii++) {
    xnames[ii] = base->xnames[ii];
  }
  for(ii=0;ii<egfs_model->np;ii++) {
    xnames[ii+base->xdim] = egfs_model->keys[egfs_model->ndf+ii];
  }
  
  cmblkl_set_names(res, xnames,err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames);
  
  return res;
}
