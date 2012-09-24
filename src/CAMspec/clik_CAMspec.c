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

void camspec_extra_init_(int*, int*,int*,int*,int*,int*, double*,double*,int*,double*,double*,double*,int*,int*,int*,int*,double*,double*);
//CAMSPEC_EXTRA_INIT(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes)

cmblkl* clik_CAMspec_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  hsize_t ndum;
  int bok;
  cmblkl *cing;
  herr_t hstat;
  char likefile[2048],szfile[2048];
  int xcase;
  double xdim;
  char *xnames_def[] = {"A_ps_100","A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib","cal0", "cal1", "cal2","xi","A_ksz"};
  
  char *xnames_1[] = {"A_ps_100","A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib","cal0", "cal1", "cal2"};
  char *xnames_0[] = {"A_ps_143", "A_cib_143", "A_sz",  
                      "cal1", "cal2"};

  char *xnames_tot[300];

  int Nspec, nX, lmax_sz;
  int *lmaxX, *lminX, *np, *npt;
  double *X,*c_inv,*tsz,*ksz,*tszXcib;
  int X_sz,c_inv_sz,sz_temp_sz;
  int beam_Nspec,num_modes_per_beam,beam_lmax,cov_dim;
  double *beam_cov_inv,*beam_modes;
  int beam_cov_inv_sz,beam_modes_sz;
  int i,j,cnt;
  //int frq[] = {100,143,217};

  camspec_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"CAMspec already initialized",*err,__LINE__,NULL);
  
  hstat = H5LTget_attribute_int( group_id, ".", "Nspec",  &Nspec);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read Nspec in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  lminX = hdf5_int_attarray(group_id,cur_lkl,"lminX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  lmaxX = hdf5_int_attarray(group_id,cur_lkl,"lmaxX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);

  np = hdf5_int_attarray(group_id,cur_lkl,"np",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  npt = hdf5_int_attarray(group_id,cur_lkl,"npt",&Nspec,err);
  forwardError(*err,__LINE__,NULL);

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
  tsz = hdf5_double_datarray(group_id,cur_lkl,"tsz",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
  sz_temp_sz = -1;
  ksz = hdf5_double_datarray(group_id,cur_lkl,"ksz",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
  sz_temp_sz = -1;
  tszXcib = hdf5_double_datarray(group_id,cur_lkl,"tszXcib",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);

  hstat = H5LTget_attribute_int( group_id, ".", "beam_Nspec",  &beam_Nspec);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read beam_Nspec in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  if (beam_Nspec==0) {
    num_modes_per_beam = 0;
    beam_lmax = 0;
    cov_dim   = 0;
    beam_cov_inv = malloc_err(sizeof(double)*1,err);
    forwardError(*err,__LINE__,NULL);
    beam_modes = malloc_err(sizeof(double)*1,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    hstat = H5LTget_attribute_int( group_id, ".", "beam_lmax",  &beam_lmax);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read beam_lmax in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
    hstat = H5LTget_attribute_int( group_id, ".", "num_modes_per_beam",  &num_modes_per_beam);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read num_modes_per_beam in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
    hstat = H5LTget_attribute_int( group_id, ".", "cov_dim",  &cov_dim);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read cov_dim in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
    beam_cov_inv_sz = -1;
    beam_cov_inv = hdf5_double_datarray(group_id,cur_lkl,"beam_cov_inv",&beam_cov_inv_sz, err);
    forwardError(*err,__LINE__,NULL);
    beam_modes_sz = -1;
    beam_modes = hdf5_double_datarray(group_id,cur_lkl,"beam_modes",&beam_modes_sz, err);
    forwardError(*err,__LINE__,NULL);
    
  }
  camspec_extra_init_(&Nspec, &nX,lminX,lmaxX,np,npt,c_inv,X,&lmax_sz, tsz,ksz,tszXcib,&beam_Nspec,&num_modes_per_beam,&beam_lmax,&cov_dim,beam_cov_inv,beam_modes);    
  
  
  //camspec_extra_getcase_(&xcase);
  xdim = 13 + beam_Nspec*num_modes_per_beam;
  
  for(i=0;i<13;i++) {
    xnames_tot[i] = xnames_def[i];
  }

  cnt = 13;
  
  for(i=0;i<beam_Nspec;i++) {
    for(j=0;j<num_modes_per_beam;j++) {
      xnames_tot[cnt] = malloc_err(sizeof(char)*50,err);
      forwardError(*err,__LINE__,NULL);
      sprintf(xnames_tot[cnt],"Bm_%d_%d",i+1,j+1);
      cnt++;
    }
  }

  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  /*if (xcase==0) {
    cmblkl_set_names(cing, xnames_0,err);
  forwardError(*err,__LINE__,NULL);  
  } else {*/
  cmblkl_set_names(cing, xnames_tot,err);
  forwardError(*err,__LINE__,NULL);
  //}
  
  free(X);
  free(c_inv);
  free(ksz);
  free(tsz);
  free(tszXcib);
  free(beam_modes);
  free(beam_cov_inv);
  for(i=13;i<xdim;i++) {
    free(xnames_tot[i]);
  }

  return cing;
}
