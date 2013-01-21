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
  gibbs_extra_lkl_(&lkl,&gb->handle,pars);
  return lkl;
}

cmblkl* clik_gibbs_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
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
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  firstchain = cldf_readint(df,"firstchain",err);
  forwardError(*err,__LINE__,NULL);
  lastchain = cldf_readint(df,"lastchain",err);
  forwardError(*err,__LINE__,NULL);
  firstsample = cldf_readint(df,"firstsample",err);
  forwardError(*err,__LINE__,NULL);
  lastsample = cldf_readint(df,"lastsample",err);
  forwardError(*err,__LINE__,NULL);
  step = cldf_readint(df,"step",err);
  forwardError(*err,__LINE__,NULL);


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

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(gb, &gibbs_lkl, 
                     &free_gibbs,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
