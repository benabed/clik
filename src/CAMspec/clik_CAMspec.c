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
  
  camspec_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"CAMspec already initialized",*err,__LINE__,NULL);
  
  memset(likefile,0,2048*sizeof(char));
  hstat = H5LTget_attribute_string( group_id, ".", "likefile",   likefile);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read likefile in %s (got %d)",*err,__LINE__,,cur_lkl,hstat);
  memset(szfile,0,2048*sizeof(char));
  hstat = H5LTget_attribute_string( group_id, ".", "szfile",   szfile);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read szfile in %s (got %d)",*err,__LINE__,,cur_lkl,hstat);
    
  camspec_extra_init_(likefile,szfile);
  
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
  
  return cing;
}
