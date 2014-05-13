#include "clik_lensing.h"
#include "lenslike/plenslike/plenslike.h"

#define _dealwitherr error *lerr,**err; if(_err==NULL) {lerr=NULL;err=&lerr;} else {err=_err;}
#define _forwardError(A,B,C) if(_err!=NULL) {forwardError(A,B,C);} else {quitOnError(A,B,stderr);}
#define _testErrorRetVA(A,B,C,D,E,F,...) if(_err!=NULL) {testErrorRetVA(A,B,C,D,E,F,__VA_ARGS__);} else {testErrorExitVA(A,B,C,D,E,__VA_ARGS__);}

int find_in_file(FILE* f, char *what, error **err) {
  char *buf;
  int i,cnt,l,j;

  cnt = 0;
  l = strlen(what);
  buf = malloc_err(sizeof(char)*(2*l+1),err);
  forwardError(*err,__LINE__,0);
  memset(buf,'\0',sizeof(char)*(2*l+2));
  j=0;
  while(1) {
    i=fread((buf+cnt),sizeof(char),1,f);
    if (i==0) {
      free(buf);
      return 0;
    }
    if (cnt==2*l-1) {
      memcpy(buf,buf+l,sizeof(char)*l);
      cnt = l-1;
    }
    cnt++;
    j++;
    buf[cnt]='\0';
    if ((cnt>=l) && (strcmp(buf+(cnt-l),what)==0)) {
      free(buf);
      return 1;
    }
  } 
}


clik_lensing_object* _clik_lensing_init(char *fpath, error **err) {
  clik_lensing_object *plid;
  FILE *f;
  cldf *df;
  int good;
  char buf[1000];
  int cnt;
  char *npath,*rpath;

  plid = malloc_err(sizeof(clik_lensing_object),err);
  forwardError(*err,__LINE__,NULL);
  plid->renorm = 1;
  plid->ren1 = 1;
  plid->has_check = 0;

  for(cnt=0;cnt<7;cnt++) {
    plid->lmax[cnt] = -1;
  }

  // test whether we have a clik file
  df = cldf_open(fpath,err);
  if (isError(*err) ) {
    // regular file
    purgeError(err);
    
    f = fopen_err(fpath,"r",err);
    forwardError(*err,__LINE__,NULL);
    
    good = find_in_file(f,"# planck lensing",err);
    forwardError(*err,__LINE__,NULL);
    testErrorRetVA(good==0,-2432,"%s is not a clik lensing file\n",*err,__LINE__,NULL,fpath);

    good = find_in_file(f,"# format: ",err);
    forwardError(*err,__LINE__,NULL);

    testErrorRetVA(good==0,-2432,"%s is not a clik lensing file (can't find format)\n",*err,__LINE__,NULL,fpath);


    memset(buf,'\0',sizeof(char)*5);
    testErrorRetVA(fread(buf,sizeof(char),4,f)!=4,-2432,"%s is not a clik lensing file (can't find format)\n",*err,__LINE__,NULL,fpath);

    buf[5] = '\0';
    npath = fpath;
    plid->type=0;
    if (strcmp(buf,"mono")==0) {
      // mono file
      fclose(f);
      plid->type=1;
    }
    if (strcmp(buf,"qecl")==0) {
      fclose(f);
      plid->type=2;
    }
    if (strcmp(buf,"full")==0) {
      fclose(f);
      plid->type=3;
    }

    testErrorRetVA(plid->type==0,-2432,"%s is not a clik lensing file (can't find format %s)\n",*err,__LINE__,NULL,fpath,buf);
  } else {
    int hk;
    plid->type = cldf_readint(df,"clik_lensing/itype",err);
    forwardError(*err,__LINE__,NULL);
    int sz=-1;
    rpath = cldf_readstr(df,"clik_lensing/filename",&sz,err);
    forwardError(*err,__LINE__,NULL);
    sprintf(buf,"%s/clik_lensing/%s",fpath,rpath);
    npath = buf;
    plid->renorm = cldf_readint(df,"clik_lensing/renorm",err);
    forwardError(*err,__LINE__,NULL);
    plid->ren1 = cldf_readint_default(df,"clik_lensing/ren1",1,err);
    forwardError(*err,__LINE__,NULL);
    hk = cldf_haskey(df,"clik_lensing/check",err);
    forwardError(*err,__LINE__,NULL);
    if (hk==1) {
      plid->check = cldf_readfloat(df,"clik_lensing/check",err);
      forwardError(*err,__LINE__,NULL);
      plid->has_check = 1;
    }
  }
    
  if (plid->type==1) {
      
    plid->plens_payload = malloc_err(sizeof(plenslike_dat_mono),err);
    forwardError(*err,__LINE__,NULL);
    
    load_plenslike_dat_mono(plid->plens_payload, npath);
    plid->lmax[0] = ((plenslike_dat_mono*) (plid->plens_payload))->lmax;
    plid->lmax[1] = ((plenslike_dat_mono*) (plid->plens_payload))->lmax;
      
  }
  if (plid->type==2) {
    plid->plens_payload = malloc_err(sizeof(plenslike_dat_qecl),err);
    forwardError(*err,__LINE__,NULL);
    good = load_plenslike_dat_qecl(plid->plens_payload, npath);
    testErrorRetVA(good!=0,-2432,"can't read %s\n",*err,__LINE__,NULL,fpath);
    plid->lmax[0] = ((plenslike_dat_qecl*) (plid->plens_payload))->lmaxphi;
    plid->lmax[1] = ((plenslike_dat_qecl*) (plid->plens_payload))->lmaxcmb;
    plid->lmax[2] = plid->lmax[1];
    plid->lmax[4] = plid->lmax[1];
  }
  if (plid->type==3) {
    plid->plens_payload = malloc_err(sizeof(plenslike_dat_full),err);
    forwardError(*err,__LINE__,NULL);
    good = load_plenslike_dat_full(plid->plens_payload, npath);
    testErrorRetVA(good!=0,-2432,"can't read %s\n",*err,__LINE__,NULL,fpath);
    plid->lmax[0] = ((plenslike_dat_full*) (plid->plens_payload))->lmaxphi;
    plid->lmax[1] = ((plenslike_dat_full*) (plid->plens_payload))->lmaxcmb;
    plid->lmax[2] = plid->lmax[1];
    plid->lmax[4] = plid->lmax[1];
  }


  return plid;
}

clik_lensing_object* clik_lensing_init(char *fpath, error **_err) {
  clik_lensing_object *plid;
  _dealwitherr;

  plid = _clik_lensing_init(fpath,err);
  _forwardError(*err,__LINE__,NULL);

  clik_lensing_selftest(plid,fpath,err);
  _forwardError(*err,__LINE__,NULL);

  return plid;
}

double clik_lensing_compute(clik_lensing_object *lclik, double *pars,error **_err) {
  double *cltt,*clte,*clee, *clphi,*clbb;
  int nextra,lmax[7];
  double lkl;
  _dealwitherr;
  
   /*nextra = int clik_lensing_get_extra_parameter_names(clikid, NULL, err);
  _forwardError(*err,__LINE__,NULL);*/
  clik_lensing_get_lmaxs(lclik,lmax,err);
  _forwardError(*err,__LINE__,-1);

  clphi = pars;
  cltt = clphi + lclik->lmax[0]+1;
  clee = cltt + lclik->lmax[1]+1;
  clbb = clee + lclik->lmax[2]+1;
  clte = clbb + lclik->lmax[3]+1;

  if (lclik->type==1) {
    plenslike_dat_mono *pmono;
    pmono = lclik->plens_payload;
    if (lclik->renorm==1) {
      lkl = calc_plenslike_mono_renorm( pmono, clphi, cltt, pmono->bl_fid);    
    } else {
      lkl = calc_plenslike_mono( pmono, clphi);
    }
  }
  if (lclik->type==2) {
    plenslike_dat_qecl *pqecl;
    pqecl = lclik->plens_payload;
    if (lclik->renorm==1) {
      lkl = calc_plenslike_qecl_renorm( pqecl, clphi, cltt, clee,clte);
    } else {
      lkl = calc_plenslike_qecl( pqecl, clphi);
    }
  }
  if (lclik->type==3) {
    plenslike_dat_full *pfull;
    pfull = lclik->plens_payload;
    if (lclik->renorm==1) {
      if (lclik->ren1==1) {
        //_DEBUGHERE_("n0 n1","");
        lkl = calc_plenslike_full_renorm_ren1( pfull, clphi, cltt, clee,clte);
      } else {
        //_DEBUGHERE_("n0","");
        lkl = calc_plenslike_full_renorm( pfull, clphi, cltt, clee,clte);  
      }
    } else {
      if (lclik->ren1==1) {
        //_DEBUGHERE_("n1","");
        lkl = calc_plenslike_full_ren1( pfull, clphi);
      } else {
        //_DEBUGHERE_("","");
        lkl = calc_plenslike_full( pfull, clphi);
      }
    }
  }

  return lkl;
}

void clik_lensing_cleanup(clik_lensing_object **plclik) {
  
  if ((*plclik)->type==1) {
    plenslike_dat_mono *pmono;
    pmono = (*plclik)->plens_payload;
    free_plenslike_dat_mono(pmono);  
  }
  if ((*plclik)->type==2) {
    plenslike_dat_qecl *pqecl;
    pqecl = (*plclik)->plens_payload;
    free_plenslike_dat_qecl(pqecl);  
  }
  if ((*plclik)->type==3) {
    plenslike_dat_full *pfull;
    pfull = (*plclik)->plens_payload;
    free_plenslike_dat_full(pfull);  
  }

  free(*plclik);
  *plclik = NULL;  
}

int clik_try_lensing(char *fpath,error **_err) {
  clik_lensing_object *plid;
  _dealwitherr;

  plid = _clik_lensing_init(fpath,err);
  if (isError(*err)) {
    //printError(stderr,*err);
    purgeError(err);
  
    return 0;
  }
  clik_lensing_cleanup(&plid);
  
  return 1;
}

int clik_lensing_get_lmax(clik_lensing_object *lclik, error **_err) {
  int lmax[7];
  _dealwitherr;

  _DEBUGHERE_("WARNING< THIS FONCTION IS DEPRECATED","");
  testErrorRet(lclik->type!=1,-12134,"this deprecated function only works with 'mono' lensing likelihood types",*err,__LINE__,-1);
  clik_lensing_get_lmaxs(lclik,lmax,err);
  forwardError(*err,__LINE__,-1);
  return lmax[0];
}

void clik_lensing_get_lmaxs(clik_lensing_object *lclik, int *lmax, error **_err) {
  int i;
  _dealwitherr;

  memset(lmax,-1,sizeof(int)*7);
  lmax[0] = lclik->lmax[0];
  for(i=0;i<7*lclik->renorm;i++) {
    lmax[i] = lclik->lmax[i];
  }
}

int clik_lensing_get_extra_parameter_names(clik_lensing_object* lclik, parnam **names, error **_err) {
  parnam *pn;
  _dealwitherr;

  if (names!=NULL) {
    pn = malloc_err(1*sizeof(parnam),err);
    _forwardError(*err,__LINE__,-1);
    *names = pn;  
  }
  return 0;
}

double* clik_lensing_clcmb_fid(clik_lensing_object* lclik, error **_err) {
  
  _dealwitherr;

  double *cltt;
  int *lmax;
  int tot,i;

  lmax = lclik->lmax;

  tot=lmax[1]+1;
  for(i=2;i<7;i++) {
    tot+= lmax[i]+1;
  }

  cltt = malloc_err(sizeof(double)*(tot),err);
  _forwardError(*err,__LINE__,NULL);

  if ((lclik)->type==1) {
    plenslike_dat_mono *pmono;
    pmono = (lclik)->plens_payload;
    memcpy(cltt,pmono->cltt_fid,sizeof(double)*(tot));    
  }
  if ((lclik)->type==2) {
    plenslike_dat_qecl *pqecl;
    pqecl = (lclik)->plens_payload;
    memcpy(cltt,pqecl->cltt_fid,sizeof(double)*(lmax[1]+1));
    memcpy(cltt+(lmax[1]+1),pqecl->clee_fid,sizeof(double)*(lmax[2]+1));
    memcpy(cltt+(lmax[1]+1+lmax[2]+1+lmax[3]+1),pqecl->clte_fid,sizeof(double)*(lmax[4]+1));
  }
  if ((lclik)->type==3) {
    plenslike_dat_full *pfull;
    pfull = (lclik)->plens_payload;
    memcpy(cltt,pfull->cltt_fid,sizeof(double)*(lmax[1]+1));
    memcpy(cltt+(lmax[1]+1),pfull->clee_fid,sizeof(double)*(lmax[2]+1));
    memcpy(cltt+(lmax[1]+1+lmax[2]+1+lmax[3]+1),pfull->clte_fid,sizeof(double)*(lmax[4]+1));
  }

  return cltt;
}

double* clik_lensing_cltt_fid(clik_lensing_object* lclik, error **_err) {
  double *cltt;
  _dealwitherr;

  cltt = clik_lensing_clcmb_fid(lclik,err);
  _forwardError(*err,__LINE__,NULL);

  return cltt;
}


double* clik_lensing_clpp_fid(clik_lensing_object* lclik, error **_err) {
  double *cltt;
  int lmax[7];
  _dealwitherr;


  clik_lensing_get_lmaxs(lclik,lmax,err);
  _forwardError(*err,__LINE__,NULL);

  cltt = malloc_err(sizeof(double)*(lmax[0]+1),err);
  _forwardError(*err,__LINE__,NULL);

  if ((lclik)->type==1) {
    plenslike_dat_mono *pmono;
    pmono = (lclik)->plens_payload;
    memcpy(cltt,pmono->clpp_fid,sizeof(double)*(lmax[0]+1));    
  }
  if ((lclik)->type==2) {
    plenslike_dat_qecl *pqecl;
    pqecl = (lclik)->plens_payload;
    memcpy(cltt,pqecl->clpp_fid,sizeof(double)*(lmax[0]+1));
  }
  if ((lclik)->type==3) {
    plenslike_dat_full *pfull;
    pfull = (lclik)->plens_payload;
    memcpy(cltt,pfull->clpp_fid,sizeof(double)*(lmax[0]+1));
  }

  return cltt;
}

void clik_lensing_selftest(clik_lensing_object *lclik, char *fpath, error **err) {
  int nextra,ndim;
  double *clX, *clt;
  double res;
  int lmax[7];
  int i,cmbdim;

  clik_lensing_get_lmaxs(lclik,lmax,err);
  forwardError(*err,__LINE__,);

  nextra = clik_lensing_get_extra_parameter_names(lclik,NULL,err);
  forwardError(*err,__LINE__,);

  ndim = nextra + lmax[0]+1;
  cmbdim = 0;
  for(i=1;i<7;i++) {
    cmbdim += lmax[i]+1;
  }
  ndim +=cmbdim;

  clt = malloc_err(sizeof(double)*ndim,err);
  forwardError(*err,__LINE__,);

  clX = clik_lensing_clpp_fid(lclik,err);
  forwardError(*err,__LINE__,);

  memcpy(clt,clX,sizeof(double)*(lmax[0]+1));

  free(clX);
  
  if (cmbdim>0) {
  clX = clik_lensing_clcmb_fid(lclik,err);
  forwardError(*err,__LINE__,);

  memcpy(clt+(lmax[0]+1),clX,sizeof(double)*(ndim-nextra-(lmax[0]+1)));

  free(clX);
  }

  res = clik_lensing_compute(lclik,clt,err);
  forwardError(*err,__LINE__,);
  if (lclik->has_check==1) {
    printf("Checking lensing likelihood '%s' on test data. got %g (expected %g)\n",fpath,res,lclik->check);  
  }
  printf("Checking lensing likelihood '%s' on test data. got %g\n",fpath,res);

  free(clt);
}

