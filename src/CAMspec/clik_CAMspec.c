#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

void free_CAMspec(void **none) {
  camspec_extra_free_();
}

double CAMspec_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  camspec_extra_lkl_(&lkl,pars);
  return lkl;
}

void camspec_extra_init_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,double*,double*,int*,double*);

cmblkl* clik_CAMspec_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  hsize_t ndum;
  int bok;
  cmblkl *cing;
  herr_t hstat;
  char likefile[2048],szfile[2048];
  int xcase;
  double xdim;
  char *xnames_1[] = {"A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib", "cal1", "cal2"};
  char *xnames_0[] = {"A_ps_143", "A_cib_143", "A_sz",  
                      "cal1", "cal2"};
  int lmin_143x143,lmax_143x143,lmin_217x217,lmax_217x217,lmin_143x217,lmax_143x217,np_143x143,np_217x217,np_143x217,nX,X_sz,c_inv_sz,sz_temp_sz,lmax_sz;
  double *X,*c_inv,*sz_temp;

  camspec_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"CAMspec already initialized",*err,__LINE__,NULL);
  
  hstat = H5LTget_attribute_int( group_id, ".", "lmin_143x143",  &lmin_143x143);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin_143x143 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax_143x143",  &lmax_143x143);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax_143x143 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "np_143x143",  &np_143x143);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read np_143x143 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "lmin_217x217",  &lmin_217x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin_217x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax_217x217",  &lmax_217x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax_217x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "np_217x217",  &np_217x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read np_217x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "lmin_143x217",  &lmin_143x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin_143x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax_143x217",  &lmax_143x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax_143x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "np_143x217",  &np_143x217);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read np_143x217 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "nX",  &nX);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read nX in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "lmax_sz",  &lmax_sz);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax_sz in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  X_sz = -1;
  X = hdf5_double_datarray(group_id,cur_lkl,"X",&X_sz, err);
  forwardError(*err,__LINE__,NULL);
  
  c_inv_sz = -1;
  c_inv = hdf5_double_datarray(group_id,cur_lkl,"c_inv",&c_inv_sz, err);
  forwardError(*err,__LINE__,NULL);

  sz_temp_sz = -1;
  sz_temp = hdf5_double_datarray(group_id,cur_lkl,"sz_temp",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
    
  camspec_extra_init_(&lmin_143x143,&lmax_143x143,&lmin_217x217,&lmax_217x217,&lmin_143x217,&lmax_143x217,&np_143x143,&np_217x217,&np_143x217,&nX,X,c_inv,&lmax_sz,sz_temp);

  
  camspec_extra_getcase_(&xcase);
  xdim = 5+4*xcase;
  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  if (xcase==0) {
    cmblkl_set_names(cing, xnames_0,err);
  forwardError(*err,__LINE__,NULL);  
  } else {
    cmblkl_set_names(cing, xnames_1,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  free(X);
  free(c_inv);
  free(sz_temp);
  return cing;
}
