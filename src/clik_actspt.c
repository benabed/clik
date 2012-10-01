#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

typedef struct {
  char tmpdir[800];
  } actspt;


void free_actspt(void **none) {
  actspt_extra_free_();
}

double actspt_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  actspt_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_actspt_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  hsize_t ndum;
  char directory_name[4096],pwd[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  herr_t hstat;
  int ilmin11,ilmin12,ilmin22,ilmax11,ilmax12,ilmax22,itt_lmax_mc,iuse_act_south  , iuse_act_equa   , iuse_spt_lowell , iuse_spt_highell;
  char dir_data[2048];
  int ldd;
  int xdim;
  char *xnames_def[] = {"a_tsz","a_ksz", "xi", "a_ps_148","a_ps_217","a_ps_95","a_ps_150","a_ps_220","a_cib_150",
        "a_cib_220","r_ps_0","r_ps_1","r_ps","r_cib","c_as_1","c_as_2","c_ae_1","c_ae_2","cal_1","cal_2","cal_3"};


  
  actspt_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"actspt already initialized",*err,__LINE__,NULL);
  
  // get data and change dir
  clik_external_data_init(directory_name,pwd,group_id,cur_lkl,err);
  forwardError(*err,__LINE__,NULL);
  
  hstat = H5LTget_attribute_int( group_id, ".", "lmin11",  &ilmin11);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin11 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmin12",  &ilmin12);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin12 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmin22",  &ilmin22);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin22 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax11",  &ilmax11);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax11 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax12",  &ilmax12);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax12 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lmax22",  &ilmax22);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax22 in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);


  hstat = H5LTget_attribute_int( group_id, ".", "tt_lmax_mc",  &itt_lmax_mc);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read tt_lmax_mc in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "use_act_equa",  &iuse_act_equa);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read use_act_equa in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "use_spt_lowell",  &iuse_spt_lowell);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read use_spt_lowell in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "use_act_south",  &iuse_act_south);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read use_act_south in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  hstat = H5LTget_attribute_int( group_id, ".", "use_spt_highell",  &iuse_spt_highell);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read use_spt_highell in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);


  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,"data/");
  dir_data[5] = ' ';
  ldd = 5;
  // call actspt_init
  actspt_extra_parameter_init_(dir_data,&ldd,&ilmin11,&ilmin12,&ilmin22,&ilmax11,&ilmax12,&ilmax22,&itt_lmax_mc,&iuse_act_south  , &iuse_act_equa   , &iuse_spt_lowell , &iuse_spt_highell);

  clik_external_data_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
    
  xdim = 21;

  cing = init_cmblkl(NULL, &actspt_lkl, 
                     &free_actspt,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);

  cmblkl_set_names(cing, xnames_def,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
