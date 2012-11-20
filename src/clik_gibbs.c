#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

typedef struct {
  int handle;
  } gibbs;


void free_gibbs(void **none) {
  gibbs *gb;

  
  gb = *none;
  gibbs_extra_free_(&(gb->handle));
  
}

double gibbs_lkl(void* none, double* pars, error **err) {
  double lkl;
  gibbs *gb;

  gb = none;
  _DEBUGHERE_("%g %g %g %g",pars[0],pars[1],pars[2],pars[3]);
  gibbs_extra_lkl_(&lkl,&gb->handle,pars);
  return lkl;
}

cmblkl* clik_gibbs_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  hsize_t ndum;
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  herr_t hstat;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step;
  gibbs *gb;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  clik_external_data_init(directory_name,pwd,group_id,cur_lkl,err);
  forwardError(*err,__LINE__,NULL);
  
  hstat = H5LTget_attribute_int( group_id, ".", "firstchain",  &firstchain);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read firstchain in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lastchain",  &lastchain);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lastchain in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "firstsample",  &firstsample);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read firstsample in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "lastsample",  &lastsample);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lastsample in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "step",  &step);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read step in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,"data/");
  dir_data[5] = ' ';
  ldd = 5;
  
  gb = malloc_err(sizeof(gibbs),err);
  forwardError(*err,__LINE__,NULL);
  
  gb->handle=0;


  //call
  gibbs_extra_parameter_init_(&(gb->handle),dir_data,&ldd,&lmin,&lmax,&firstchain,&lastchain,&firstsample,&lastsample,&step);
  testErrorRetVA(gb->handle<=0,-4325325432,"handle return is negative (got %d)",*err,__LINE__,NULL,gb->handle);

  clik_external_data_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(gb, &gibbs_lkl, 
                     &free_gibbs,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
