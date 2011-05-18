#include "clik.h"
#include "hdf5.h"
#include "lklbs.h"
#include "aplowly.h"
#include "fowly.h"
#include "smica.h"
#include <errno.h>
#include <string.h>

typedef struct {
  char tmpdir[800];
  } bopix;


void free_bopix(void **none) {
  bopix_extra_free_();
}

double bopix_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  bopix_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_bopix_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char dirtmpl[2048];
  char *drn;
  hsize_t ndum;
  char *bopix_data;
  char fpix_data_name[2048];
  FILE *fpix_data;
  char command[4096];
  int status;
  char *bopix_nthreads_env;
  int bopix_nthreads;
  int bok;
  cmblkl *cing;
  int mlmax;
  herr_t hstat;
  
  bopix_extra_only_one_(&bok);
  testErrorRet(bok==0,-100,"Bopix already initialized",*err,__LINE__,NULL);
  
  // create temporary dir
  sprintf(dirtmpl,"/tmp/bopix_XXXXXX");
  drn = mkdtemp(dirtmpl);
  testErrorRetVA(drn==NULL,-100,"cannot create temporary dir (cause = '%s')",*err,__LINE__,NULL,strerror(errno));
  
  // read tarfile from hdffile
  hstat = H5LTget_dataset_info( group_id, "bopix_data", &ndum, NULL,NULL);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,"bopix_data",cur_lkl,hstat);
  bopix_data = malloc_err(sizeof(char)*ndum,err);
  forwardError(*err,__LINE__,NULL);
  hstat = H5LTread_dataset_double(group_id,"bopix_data",bopix_data);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,"bopix_data",cur_lkl,hstat);
  
  // save to file !
  sprintf(fpix_data_name,"%s/bopix_data.tar",drn);
  fpix_data = fopen_err(fpix_data_name,"w",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(fwrite(bopix_data,1,ndum,fpix_data)<ndum,-100,"Cannot write to file %s",*err,__LINE__,NULL,fpix_data_name);
  fclose(fpix_data);
  free(bopix_data);
  
  // change dir
  testErrorRetVA(chdir(drn)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,NULL,drn,strerror(errno));

  // call tar to recreate the files  
  sprintf(command,"tar xf %s",fpix_data_name);
  status = system(command);
  testErrorRetVA(status!=0,-100,"cannot untar, command '%s' got status %d",*err,__LINE__,NULL,command,status);
  
  // call bopix_init
  bopix_extra_parameter_init_();
  
  // hope for the best !
  
  // delete all files (like a macho !)
  sprintf(command,"rm -rf %s",drn); 
  status = system(command);
  testErrorRetVA(status!=0,-100,"cannot delete files, command '%s' got status %d",*err,__LINE__,NULL,command,status);
  
  bopix_nthreads = clik_getenviron_numthread("BOPIX",1,err);
  forwardError(*err,__LINE__,NULL);
  bopix_extra_set_bopix_nthreads_(&bopix_nthreads);
  
  mlmax = ell[nell-1];
  bopix_extra_set_lmax_(&mlmax);

  cing = init_cmblkl(NULL, &bopix_lkl, 
                     &free_bopix,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;
}
