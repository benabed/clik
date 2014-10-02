#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

typedef struct {
  char tmpdir[800];
  } plik_cmbonly;


void free_plik_cmbonly(void **none) {
  plik_cmbonly_extra_free_();
}

double plik_cmbonly_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  //_DEBUGHERE_("%g %g %g %g",pars[0],pars[1],pars[2],pars[3]);
  plik_cmbonly_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_plik_cmbonly_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  

  plik_cmbonly_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"plik_cmbonly already initialized",*err,__LINE__,NULL);
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
    
  memset(dir_data,' ',sizeof(char)*2048);
  //sprintf(dir_data,"");
  //dir_data[5] = ' ';
  ldd = 0;
  
  // call plik_cmbonly_init
  plik_cmbonly_extra_init_(dir_data,&ldd);

  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(NULL, &plik_cmbonly_lkl, 
                     &free_plik_cmbonly,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);


  return cing;
}
