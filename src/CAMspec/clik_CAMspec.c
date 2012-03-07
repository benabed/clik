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

void camspec_extra_init_(int*, int*,int*,int*,int*,int*,int*, double*,double*,double*,double*);
//CAMSPEC_EXTRA_INIT(iNspec, inX,ilminX,ilmaxX,inp,inpt,ilmax_sz, sz_100,sz_143,mc_inv,mX)

cmblkl* clik_CAMspec_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  hsize_t ndum;
  int bok;
  cmblkl *cing;
  herr_t hstat;
  char likefile[2048],szfile[2048];
  int xcase;
  double xdim;
  char *xnames_1[] = {"A_ps_100","A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib","cal0", "cal1", "cal2"};
  char *xnames_0[] = {"A_ps_143", "A_cib_143", "A_sz",  
                      "cal1", "cal2"};
  int Nspec, nX, lmax_sz;
  int *lmaxX, *lminX, *np, *npt;
  double *X,*c_inv,*sz_143,*sz_100;
  int X_sz,c_inv_sz,sz_temp_sz;

  _DEBUGHERE_("","");
  camspec_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"CAMspec already initialized",*err,__LINE__,NULL);
  
  _DEBUGHERE_("","");
  hstat = H5LTget_attribute_int( group_id, ".", "Nspec",  &Nspec);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read Nspec in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  _DEBUGHERE_("","");
  lminX = hdf5_int_attarray(group_id,cur_lkl,"lminX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  lmaxX = hdf5_int_attarray(group_id,cur_lkl,"lmaxX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);

  _DEBUGHERE_("","");
  np = hdf5_int_attarray(group_id,cur_lkl,"np",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  npt = hdf5_int_attarray(group_id,cur_lkl,"npt",&Nspec,err);
  forwardError(*err,__LINE__,NULL);

  _DEBUGHERE_("","");
  hstat = H5LTget_attribute_int( group_id, ".", "nX",  &nX);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read nX in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  _DEBUGHERE_("","");
  hstat = H5LTget_attribute_int( group_id, ".", "lmax_sz",  &lmax_sz);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax_sz in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  _DEBUGHERE_("","");
  X_sz = -1;
  X = hdf5_double_datarray(group_id,cur_lkl,"X",&X_sz, err);
  forwardError(*err,__LINE__,NULL);
  
  _DEBUGHERE_("","");
  c_inv_sz = -1;
  c_inv = hdf5_double_datarray(group_id,cur_lkl,"c_inv",&c_inv_sz, err);
  forwardError(*err,__LINE__,NULL);

  _DEBUGHERE_("","");
  sz_temp_sz = -1;
  sz_143 = hdf5_double_datarray(group_id,cur_lkl,"sz_143",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
  
  _DEBUGHERE_("","");
  sz_temp_sz = -1;
  sz_100 = hdf5_double_datarray(group_id,cur_lkl,"sz_100",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);

  _DEBUGHERE_("","");
  _DEBUGHERE_("%d %d %d %d",lminX[0],lminX[1],lminX[2],lminX[3])
  camspec_extra_init_(&Nspec, &nX,lminX,lmaxX,np,npt,&lmax_sz, sz_100,sz_143,c_inv,X);    
  
  _DEBUGHERE_("","");
  
  //camspec_extra_getcase_(&xcase);
  xdim = 11;
  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  /*if (xcase==0) {
    cmblkl_set_names(cing, xnames_0,err);
  forwardError(*err,__LINE__,NULL);  
  } else {*/
  _DEBUGHERE_("","");
  cmblkl_set_names(cing, xnames_1,err);
  forwardError(*err,__LINE__,NULL);
  //}
  _DEBUGHERE_("","");
  
  free(X);
  free(c_inv);
  free(sz_143);
  free(sz_100);
  _DEBUGHERE_("","");
  return cing;
}
