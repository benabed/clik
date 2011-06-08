#include "clik_helper.h"

// NO USER SERVICEABLE PART HERE.
// THINGS CAN CHANGE WITHOUT NOTICE. 



// ARE YOU STILL READING ?

// YOU HAVE BEEN WARNED !


#define _dealwitherr error *lerr,**err; if(_err==NULL) {lerr=NULL;err=&lerr;} else {err=_err;}

#define _forwardError(A,B,C) if(_err!=NULL) {forwardError(A,B,C);} else {quitOnError(A,B,stderr);}
#define _testErrorRetVA(A,B,C,D,E,F,...) if(_err!=NULL) {testErrorRetVA(A,B,C,D,E,F,__VA_ARGS__);} else {testErrorExitVA(A,B,C,D,E,__VA_ARGS__);}

double* hdf5_double_datarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err) {
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  double *res;
  
  hstat = H5LTget_dataset_info( group_id, name, &ndum, &dum, &ddum);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  testErrorRetVA((ndum!=*sz && *sz>0),hdf5_base,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,name,cur_lkl,ndum,*sz);
  res = malloc_err(sizeof(double)*ndum,err);
  forwardError(*err,__LINE__,NULL);
  hstat = H5LTread_dataset_double(group_id, name,res);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  if (*sz<0) {
    *sz = ndum;
  }
  return res;
}

int* hdf5_int_datarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err) {
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  int *res;
  
  hstat = H5LTget_dataset_info( group_id, name, &ndum, &dum, &ddum);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  testErrorRetVA((ndum!=*sz && *sz>0),hdf5_base,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,name,cur_lkl,ndum,*sz);
  res = malloc_err(sizeof(int)*ndum,err);
  forwardError(*err,__LINE__,NULL);
  hstat = H5LTread_dataset_int(group_id, name,res);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  if (*sz<0) {
    *sz = ndum;
  }
  return res;
  //_DEBUGHERE_("%d -> %d",res[0],res[ndum-1]);
}

double* hdf5_double_attarray(hid_t group_id,char*  cur_lkl,char* name,int* sz, error **err) {
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  double *res;
  
  hstat = H5LTget_attribute_info( group_id, ".",name, &ndum, &dum, &ddum);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  testErrorRetVA((ndum!=*sz && *sz>0),hdf5_base,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,name,cur_lkl,ndum,*sz);
  res = malloc_err(sizeof(double)*ndum,err);
  forwardError(*err,__LINE__,NULL);
  hstat = H5LTget_attribute_double(group_id, ".",name,res);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,name,cur_lkl,hstat);
  if (*sz<0) {
    *sz = ndum;
  }
  return res;
}

int clik_getenviron_integer(char* name, int sfg, error **err) {
  int res;
  char *cres;
  int flg;
  
  cres = getenv(name);
  if (cres!=NULL) {
    flg = sscanf(cres,"%d",&res);
    if (flg==1) {
      return res;
    }
  }
  return sfg;
}

double clik_getenviron_real(char* name, double sfg, error **err) {
  double res;
  char *cres;
  int flg;

  cres = getenv(name);
  if (cres!=NULL) {
    flg = sscanf(cres,"%lg",&res);
    if (flg==1) {
      return res;
    }
  }
  return sfg;
}


char* clik_getenviron_string(char* name, char* sfg, error **err) {
  char *cres;
  
  cres = getenv(name);
  if (cres!=NULL) {
    return cres;
  }
  return sfg;
}

int clik_getenviron_numthread(char* name, int sfg, error **err) {
  int np;
  char fullname[2048];
  int i;
  
  sprintf(fullname,"%s_NUMTHREADS",name);
  for(i=0;i<strlen(name);i++) {
    fullname[i] = toupper(fullname[i]);
  }
  
  np = clik_getenviron_integer(fullname,sfg,err);
  forwardError(*err,__LINE__,sfg);
  
  testErrorRetVA(np<0 && np!=sfg, -100,"%s env variable meaningless",*err,__LINE__,sfg,fullname);
  return np;
}

cmblkl * clik_lklobject_init(hid_t group_id,char* cur_lkl,error **err) {
  cmblkl *clkl;
  parname lkl_type;
  herr_t hstat;
  int has_cl[6];
  int nell, *ell,nbins,i,cli;
  double *wl,*bins;
  double unit;
  int lmin, lmax;
  char init_func_name[2048];
  clik_lkl_init_func *clik_dl_init;
  void* dlhandle;   
  
  // get the lkl type
  memset(lkl_type,0,_pn_size*sizeof(char));

  hstat = H5LTget_attribute_string( group_id, ".", "lkl_type",  lkl_type);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read lkl_type in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  // get unit
  hstat = H5LTget_attribute_double( group_id, ".", "unit",  &unit);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read unit in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  // get has_cl
  hstat = H5LTget_attribute_int( group_id, ".", "has_cl",  has_cl);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read has_cl in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  // get ells
  hstat = H5LTfind_attribute (group_id, "lmax");
  if (hstat==1) {
    // has lmax !
    hstat = H5LTget_attribute_int( group_id, ".", "lmax", &lmax);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read lmax in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
    
    lmin = 0;
    hstat = H5LTfind_attribute (group_id, "lmin");
    if (hstat==1) {
      hstat = H5LTget_attribute_int( group_id, ".", "lmin", &lmin);
      testErrorRetVA(hstat<0,hdf5_base,"cannot read lmin in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
      
    }
    nell = lmax+1-lmin;
    ell = malloc_err(sizeof(int)*(nell),err);
    forwardError(*err,__LINE__,NULL);
    for(i=lmin;i<lmax+1;i++) {
      ell[i-lmin] = i;
    }
  } else {
    hstat = H5LTget_attribute_info(group_id, ".", "ell",&nell, NULL, NULL);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read ell in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
    ell = malloc_err(sizeof(int)*(nell),err);
    forwardError(*err,__LINE__,NULL);
    hstat = H5LTget_attribute_int(group_id, ".", "ell",ell);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read ell in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  }  
  
  lmax = ell[nell-1];
  
  // get wl
  wl = NULL;
  hstat = H5LTfind_attribute (group_id, "wl");
  if (hstat==1) {
    int nwl;
    nwl = lmax+1;
    wl = hdf5_double_attarray(group_id,cur_lkl,"wl",&nwl,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  // deals with bins
  nbins = 0;
  bins = NULL;
  hstat = H5LTfind_attribute(group_id, "nbins");
  
  if (hstat==1) {
    hstat = H5LTget_attribute_int(group_id, ".", "nbins",&nbins);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read nbins in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  }
  if (nbins!=0) {
    int nd;
    int ncl;
    hsize_t nn;
    H5T_class_t dum;
    size_t ddum;
    hid_t dataset_id;
    
    nd = 0;
    ncl = 0;
    for(cli=0;cli<6;cli++) {
      if (has_cl[cli]==1) {
        nd += nell;
        ncl++;
      }
    }
    hstat = H5LTfind_dataset(group_id, "nbins");
    
    if (hstat==1) { // full binning matrix
      int nbn;
      nbn = nbins*nd;
      
      bins = hdf5_double_datarray(group_id, cur_lkl,"bins",&nbn,err);
      forwardError(*err,__LINE__,NULL);  
      
      
    } else { //packed binning matrix
      double *bin_ws;
      int *ellmin,*ellmax;
      int nw; 
      int ib,jb,il;
      int wsz;
      
      nw=-1;
      bin_ws = hdf5_double_datarray(group_id, cur_lkl,"bin_ws",&nw,err);
      forwardError(*err,__LINE__,NULL);  
      
      ellmin = hdf5_int_datarray(group_id, cur_lkl,"bin_lmin",&nbins,err);
      forwardError(*err,__LINE__,NULL);  
      ellmax = hdf5_int_datarray(group_id, cur_lkl,"bin_lmax",&nbins,err);
      forwardError(*err,__LINE__,NULL);  
      
      bins = malloc_err(sizeof(double)*nd*nbins,err);
      forwardError(*err,__LINE__,NULL);  
      memset(bins,0,sizeof(double)*nd*nbins);
      
      jb=0;
      for(ib=0;ib<nbins;ib++) {        
        wsz = ellmax[ib] - ellmin[ib] + 1;
        testErrorRetVA(jb+wsz>nw,-11111,"argl bins",*err,__LINE__,NULL,"");
        memcpy(&(bins[ib*nd+ellmin[ib]]),&(bin_ws[jb]),wsz*sizeof(double));
        jb+=wsz;
      }
      free(ellmin);
      free(ellmax);
      free(bin_ws);
    }
  }
  
  clkl = NULL;

  sprintf(init_func_name,"clik_%s_init",lkl_type);
#ifdef HAS_RTLD_DEFAULT 
  dlhandle = RTLD_DEFAULT;
#else
  dlhandle = NULL;
#endif
  clik_dl_init = dlsym(dlhandle,init_func_name);
  testErrorRetVA(clik_dl_init==NULL,-1111,"Cannot initialize lkl type %s from %s dl error : %s",*err,__LINE__,NULL,lkl_type,cur_lkl,dlerror());  \

  clkl = clik_dl_init(group_id,cur_lkl,nell,ell,has_cl,unit,wl,bins,nbins,err);
  forwardError(*err,__LINE__,NULL); 

  /*if (strcasecmp(lkl_type, "ivg")==0) {
    clkl = clik_ivg_init(group_id,cur_lkl,nell,ell,has_cl,unit,wl,bins,nbins,err);
    forwardError(*err,__LINE__,NULL);
  }  
  if (strcasecmp(lkl_type, "gauss")==0) {
    clkl = clik_gauss_init(group_id,cur_lkl,nell,ell,has_cl,unit,wl,bins,nbins,err);
    forwardError(*err,__LINE__,NULL);
  }
  if (strcasecmp(lkl_type, "smica")==0) {
    clkl = clik_smica_init(group_id,cur_lkl,nell,ell,has_cl,unit,wl,bins,nbins,err);
    forwardError(*err,__LINE__,NULL);
  }
  testErrorRetVA(clkl == NULL,-1000,"Unknown lkl %s",*err,__LINE__,NULL,lkl_type);
  */
  // cleanups
  if(wl!=NULL) {
    free(wl);    
  }
  if(bins!=NULL) {
    free(bins);    
  }
  free(ell);
  
  return clkl;
}

void clik_external_data_init(char *pwd,char *dirname,hid_t group_id, char* cur_lkl,error **err) {
  herr_t hstat;
  char dirtmpl[2048*4];
  char *drn;
  hsize_t ndum;
  char fpix_data_name[2048*4];
  FILE *fpix_data;
  char command[4096*4];
  int status;
   
  testErrorRetVA(getcwd(pwd,4096)==NULL,-101010,"can't get cwd name (cause = '%s')",*err,__LINE__,,strerror(errno));
  
  // do we need to extract the data ?
  hstat = H5LTfind_dataset (group_id, "external_data");
  if (hstat==1) {
    char *data;
    
    // yes !
    sprintf(dirtmpl,"/tmp/clik_XXXXXX");
    drn = mkdtemp(dirtmpl);
    testErrorRetVA(drn==NULL,-100,"cannot create temporary dir (cause = '%s')",*err,__LINE__,,strerror(errno));

    // read tarfile from hdffile
    hstat = H5LTget_dataset_info( group_id, "external_data", &ndum, NULL,NULL);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,,"tardata",cur_lkl,hstat);
    data = malloc_err(sizeof(char)*ndum,err);
    forwardError(*err,__LINE__,);
    hstat = H5LTread_dataset_char(group_id,"external_data",data);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,,"tardata",cur_lkl,hstat);

    // save to file !
    sprintf(fpix_data_name,"%s/data.tar",drn);
    fpix_data = fopen_err(fpix_data_name,"w",err);
    forwardError(*err,__LINE__,);
    testErrorRetVA(fwrite(data,1,ndum,fpix_data)<ndum,-100,"Cannot write to file %s",*err,__LINE__,,fpix_data_name);
    fclose(fpix_data);
    free(data);

    // change dir
    testErrorRetVA(chdir(drn)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,drn,strerror(errno));

    // call tar to recreate the files  
    sprintf(command,"tar xf %s",fpix_data_name);
    status = system(command);
    testErrorRetVA(status!=0,-100,"cannot untar, command '%s' got status %d",*err,__LINE__,,command,status);
    sprintf(dirname,"%s",drn);
    
  }  else {
    memset( fpix_data_name,0,2048*4);
    hstat = H5LTget_attribute_string( group_id, ".", "external_dir",   fpix_data_name);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read external_dir in %s (got %d)",*err,__LINE__,,cur_lkl,hstat);
    testErrorRetVA(chdir(fpix_data_name)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,fpix_data_name,strerror(errno));
    dirname[0]='\0';
  }

}

void clik_external_data_cleanup(char* pwd,char *dirname,error **err) {
  char command[4096*4];
  int status;
  
   // delete all files (like a macho !)
  testErrorRetVA(chdir(pwd)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,pwd,strerror(errno));
  
  if (dirname[0]!='\0') {
    // remove files
    sprintf(command,"rm -rf %s",dirname); 
    status = system(command);
    testErrorRetVA(status!=0,-100,"cannot delete files, command '%s' got status %d",*err,__LINE__,,command,status);    
  }
}