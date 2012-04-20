#include "clik_parametric.h"
#include "clik_helper.h"

typedef struct {
  parametric *p_model;
  double *rq,*wrq;
  int nell,nbins,m;
  double unit;
  int *ell;
  double *bins,*wl;
  double *A;
  } parametric_smica;

void comp_parametric_update(void* data,double* locpars, double* rq, error **err) {
  parametric_smica *p_pay;
  SmicaComp* SC;
  double *wl,*wl0,one;
  int inc,il,im;
  double res;
  //double r10[16];

  SC = data;
  p_pay = SC->data;

  //_DEBUGHERE_("%g",locpars[0]);

  parametric_compute(p_pay->p_model, locpars, p_pay->rq, NULL, err);
  forwardError(*err,__LINE__,);
  
  
  // apply wl and binning
  one=1;
  if (p_pay->wl==NULL) {
    wl0 = &one;
    inc = 0;
  } else {
    wl0 = p_pay->wl;
    inc = 1;
  }
  
  //_DEBUGHERE_("","");
  wl=wl0;
  for(il=0;il<p_pay->nell;il++) {
    int ip;
    int im1,im2;
    ip = il*p_pay->m*p_pay->m;

    for(im1=0;im1<p_pay->m;im1++) {
      for(im2=0;im2<p_pay->m;im2++) {
        p_pay->rq[il*p_pay->m*p_pay->m+im1*p_pay->m+im2] = p_pay->rq[il*p_pay->m*p_pay->m+im1*p_pay->m+im2] * *wl * p_pay->unit * p_pay->A[im1]*p_pay->A[im2];  
      }
    }
    wl+=inc;
  }
  //_DEBUGHERE_("%g %g %g %g",egfs_pay->A[0],egfs_pay->A[1],egfs_pay->A[2],egfs_pay->A[3])
  //_DEBUGHERE_("%g %g %g",egfs_pay->rq[0],egfs_pay->rq[2],egfs_pay->unit);
  
  // apply binning if needed
  if (p_pay->bins!=NULL) {
    char transa,transb;
    int npar;
    double done,dzero;
    int one,nbns,nell;
    int ndim;
  //  int ii;
  //  double rq0,rq2;
  //  double poc;

    transa='N';
    transb='N';
    ndim = p_pay->m*p_pay->m;
    one = 1;
    done = 1.;
    dzero = 0;
    nbns = p_pay->nbins;
    nell = p_pay->nell;
    
    /*rq0=rq[0];
    rq2=rq[2];
    _DEBUGHERE_("%g %g",rq[0],rq[2]);
    _DEBUGHERE_("%g %g",rq[0]-rq0,rq[2]-rq2);
    _DEBUGHERE_("m %d n %d k %d",ndim, nbns,nell);*/  
    //_DEBUGHERE_("avant egfs","");
    //memset(r10,0,sizeof(double)*16);

    //printMat(&rq[10*ndim], egfs_pay->m,egfs_pay->m);
    //printMat(r10, egfs_pay->m,egfs_pay->m);
    {
      int il,iq,if1,if2;
      for(il=0;il<nell;il++) {
        for(iq=0;iq<nbns;iq++) {
          for(if1=0;if1<p_pay->m;if1++) {
            for(if2=0;if2<p_pay->m;if2++) {
              rq[iq*ndim+if1*p_pay->m+if2] += p_pay->bins[iq*nell+il] * p_pay->rq[il*p_pay->m*p_pay->m+if1*p_pay->m+if2];
              /*if (iq==10)
                r10[if1*egfs_pay->m+if2] += egfs_pay->bins[iq*nell+il] * egfs_pay->rq[il+if1*egfs_pay->m*egfs_pay->nell+if2*egfs_pay->nell];*/
            }  
          }
        }
      }
    }
    //dgemm(&transa, &transb, &ndim, &nbns,&nell, &done, egfs_pay->rq, &ndim, egfs_pay->bins, &nell, &done, rq, &ndim);
    /*_DEBUGHERE_("apres egfs","");
    printMat(&rq[10*ndim], egfs_pay->m,egfs_pay->m);
    printMat(r10, egfs_pay->m,egfs_pay->m);*/
    /*_DEBUGHERE_("","");
    poc = 0;
    for(ii=0;ii<10;ii++) {
      _DEBUGHERE_("%g %g %g",egfs_pay->rq[ii*ndim],egfs_pay->rq[ii*ndim+2],egfs_pay->bins[ii]);
      poc += egfs_pay->rq[ii*ndim] * egfs_pay->bins[ii];
    }
    _DEBUGHERE_("%g %g %g",rq[0]-rq0,rq[2]-rq2,poc);*/
  } else {
    int if1,if2;
    for(il=0;il<p_pay->nell;il++) {
      for(if1=0;if1<p_pay->m;if1++) {
        for(if2=0;if2<p_pay->m;if2++) {
          rq[il*p_pay->m*p_pay->m+if1*p_pay->m+if2] += p_pay->rq[il*p_pay->m*p_pay->m+if1*p_pay->m+if2];
        }
      }
    }
  }
    
}
void free_comp_parametric(void** data) {
  SmicaComp *SC;
  parametric_smica *p_pay;
  
  SC = *data;
  p_pay = SC->data;
  free(p_pay->rq);
  if (p_pay->nbins!=0) {
    free(p_pay->bins);
  }
  if (p_pay->wl!=NULL) {
    free(p_pay->wl);
  }
  parametric_free(&(p_pay->p_model));
  free(p_pay->A);

  free(SC->data);
  free(SC);
  *data = NULL;
}

void base_parametric_hdf5_init(hid_t comp_id,char* cur_lkl,int ndet, double** detlist,int *ndef, char ***defkeys, char*** defvalues, int *nvar, char ***varkeys, error **err) {
  int dz,i;
  char *keyvartable,*deftable,*valtable;
  herr_t hstat;

  hstat = H5LTget_attribute_int( comp_id, ".", "ndim",  nvar);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read ndim in %s (got %d)",*err,__LINE__,,cur_lkl,hstat);

  dz = -1;
  keyvartable = hdf5_char_attarray(comp_id,cur_lkl,"keys",&dz, err);
  forwardError(*err,__LINE__,);
  
  if (*nvar!=0) {
    *varkeys = malloc_err(sizeof(char*)**nvar,err);
    forwardError(*err,__LINE__,);
  } else {
    *varkeys = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,);
    (*varkeys)[0] = NULL;
  }

  for(i=0;i<*nvar;i++) {
    (*varkeys)[i] = &(keyvartable[i*256]);
  }

  // get defaults
  hstat = H5LTget_attribute_int( comp_id, ".", "ndef",  ndef);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read ndef in %s (got %d)",*err,__LINE__,,cur_lkl,hstat);
  
  dz = -1;
  deftable = hdf5_char_attarray(comp_id,cur_lkl,"defaults",&dz, err);
  forwardError(*err,__LINE__,);
  dz = -1;
  valtable = hdf5_char_attarray(comp_id,cur_lkl,"values",&dz, err);
  forwardError(*err,__LINE__,);
  
  if (*ndef!=0) {
    *defkeys = malloc_err(sizeof(char*)**ndef,err);
    forwardError(*err,__LINE__,);
    *defvalues = malloc_err(sizeof(char*)**ndef,err);
    forwardError(*err,__LINE__,);    
  } else {
    *defkeys = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,);
    *defvalues = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,);
    (*defkeys)[0] = NULL;
    (*defvalues)[0] = NULL;
  }
  
  for(i=0;i<*ndef;i++) {
    (*defkeys)[i] = &(deftable[i*256]);
    (*defvalues)[i] = &(valtable[i*256]);
  }  
  
  dz = ndet;
  hstat = H5LTfind_attribute (comp_id, "dfreq");
  if (hstat ==1) {
    *detlist = hdf5_double_attarray(comp_id,cur_lkl,"dfreq",&dz, err);
    forwardError(*err,__LINE__,);
  } else {
    int * ietlist;
    ietlist = hdf5_int_attarray(comp_id,cur_lkl,"freq",&dz, err);
    forwardError(*err,__LINE__,);
    *detlist = malloc_err(sizeof(double)*ndet,err);
    forwardError(*err,__LINE__,);
    
    for(i=0;i<ndet;i++) {
      (*detlist)[i]=ietlist[i];
    }
    free(ietlist);
  }
  
  return;

}

SmicaComp * finalize_parametric_hdf5_init(parametric* p_model,hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  parametric_smica *p_pay;
  int i,eb;
  char **xnames;
  SmicaComp *SC;
  int lmin,lmax;
  herr_t hstat;
  double *color;
  int dz;

  hstat = H5LTfind_dataset(comp_id, "color");
  if (hstat ==1) {
    dz  = p_model->ndet;
    color = hdf5_double_datarray(comp_id,cur_lkl,"color",&dz, err);
    forwardError(*err,__LINE__,NULL);
    parametric_set_color(p_model,color,err);
    forwardError(*err,__LINE__,NULL);
    free(color);
  }
    
  lmin = ell[0];

  lmax = ell[nell-1];
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);

  eb = 0;
  for(i=1;i<6;i++) {
    eb +=has_cl[i];
  }
  testErrorRet(eb!=0,-7693,"parametric does not work with polarized data yet",*err,__LINE__,NULL);
  
  p_pay = malloc_err(sizeof(parametric_smica),err);
  forwardError(*err,__LINE__,NULL);
    
  p_pay->m = m;
  p_pay->p_model = p_model;  
  p_pay->unit = unit;
  
  p_pay->A = hdf5_double_attarray(comp_id,cur_lkl,"A_cmb",&m,err);
  forwardError(*err,__LINE__,NULL);    

  p_pay->nell = nell;

  p_pay->nbins = nbins;
  p_pay->bins = NULL;
  if (nbins !=0) {
    p_pay->bins = malloc_err(sizeof(double)*(nell*nbins),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(p_pay->bins,bins,sizeof(double)*nbins*nell);    
  }
  p_pay->wl = NULL;
  if (wl!=NULL) {
    p_pay->wl = malloc_err(sizeof(double)*(nell),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(p_pay->wl,wl,sizeof(double)*nell);    
    
  }
  p_pay->rq = malloc_err(sizeof(double)*(lmax+1-lmin)*m*m,err);
  forwardError(*err,__LINE__,NULL);
  
  SC = alloc_SC(p_model->nvar,nb,m,p_pay,&comp_parametric_update,&free_comp_parametric,err);
  forwardError(*err,__LINE__,NULL);
  
  if (p_model->nvar!=0) {
    xnames = malloc_err(sizeof(char*)*(p_model->nvar),err);
    forwardError(*err,__LINE__,NULL);
  } else{
    xnames = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,NULL);
  }   
  for(i=0;i<p_model->nvar;i++) {
    xnames[i] = p_model->varkey[i];
  }
  
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames);
  
  return SC;    
}
