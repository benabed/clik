#include "clik_helper.h"

cmblkl* clik_lowly_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  powly *self;
  cmblkl *cing;
  herr_t hstat;
  int neff,nnn,cli;
  
  self = malloc_err(sizeof(powly), err);
  forwardError(*err,__LINE__,NULL);
  SET_PRINT_STAT(self);
  
  
  hstat = H5LTget_attribute_int( group_id, ".", "neff",  &neff);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read neff in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  self->neff = neff;
  
  self->nell = nell;
  self->ell=malloc_err(sizeof(int)*nell,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(self->ell,ell,sizeof(int)*nell);
  self->ncl_fid = lowly_get_offset_cl(has_cl,self->offset_cl,self->nell);  
  self->nmode = lowly_nmode(self->nell,self->ell);
  
  self->offset_mode[0]=(-1,-1,-1);
  self->tot_mode = 0;
  if ((self->offset_cl[0]!=-1) || (self->offset_cl[3]!=-1) || (self->offset_cl[4]!=-1)) {
    self->offset_mode[0]=self->tot_mode; //T
    self->tot_mode += self->nmode;
  }
  if ((self->offset_cl[1]!=-1) || (self->offset_cl[3]!=-1) || (self->offset_cl[5]!=-1)) {
    self->offset_mode[1]=self->tot_mode; //E
    self->tot_mode += self->nmode;
  }
  if ((self->offset_cl[2]!=-1) || (self->offset_cl[4]!=-1) || (self->offset_cl[5]!=-1)) {
    self->offset_mode[1]=self->tot_mode; //B
    self->tot_mode += self->nmode;
  }
  
  self->cl_fid = hdf5_double_datarray(group_id, cur_lkl,"cl_fid",&self->ncl_fid,err);
  forwardError(*err,__LINE__,NULL);

  self->buffer = malloc_err(sizeof(double)*(_SZT_(self->tot_mode)*_SZT_(self->neff) + _SZT_(self->neff)*_SZT_(self->neff) + self->neff*2),err);
  forwardError(*err,__LINE__,NULL);
  self->H = self->buffer;
  self->a_bar = self->H + self->tot_mode*self->neff;
  self->tmp = self->a_bar + self->neff;
  self->Re = self->tmp + self->neff;
  
  nnn = self->tot_mode*self->neff;
  self->H = hdf5_double_datarray(group_id, cur_lkl,"H",&nnn,err);
  forwardError(*err,__LINE__,NULL);
  
  nnn = self->neff;
  self->a_bar = hdf5_double_datarray(group_id, cur_lkl,"a_bar",&nnn,err);
  forwardError(*err,__LINE__,NULL);
  
  self->time_build = 0;
  self->time_tot = 0;
  self->time_chol = 0;
  self->n_tests = 0;
  
  self->eMat = NULL;
  self->eHt = NULL;
  
  self->n_sym = 0;
  for(cli=0;cli<3;cli++) {
    self->n_sym += self->offset_cl[cli]==-1 ? 0:1;
  }
  self->n_cross = 0;
  for(cli=0;cli<3;cli++) {
    self->n_cross += self->offset_cl[cli+3]==-1 ? 0:1;
  }
  self ->nstore[0]=0;
  self ->nstore[1]=0;
  self ->tot_store=0;
  self->storage = NULL;
  
  cing = init_cmblkl(self, &powly_lkl, 
                     &free_powly,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}

void al_bl_cor(hid_t group_id, char* cur_lkl, int *_n, double **_al, double **_bl, double **_nl, double **_cor,int *_isdiag,error **err) {
  double *al, *bl, *corr,*nl;
  int n,isdiag;
  herr_t hstat;
  
  n = *_n;
  // get al
  al = hdf5_double_datarray(group_id, cur_lkl,"al",&n,err);
  forwardError(*err,__LINE__,);

  // get bl
  bl = hdf5_double_datarray(group_id, cur_lkl,"bl",&n,err);
  forwardError(*err,__LINE__,);

  // get nl
  nl = NULL;
  hstat = H5LTfind_dataset(group_id, "nl");
  if (hstat==1) {
    nl = hdf5_double_datarray(group_id, cur_lkl,"nl",&n,err);
    forwardError(*err,__LINE__,);
  }

  // get cor
  corr = NULL;
  hstat = H5LTfind_dataset(group_id, "cor");
  if (hstat==1) {
    int ncor;
    ncor=-1;
    corr = hdf5_double_datarray(group_id, cur_lkl,"cor",&ncor,err);
    forwardError(*err,__LINE__,);
    testErrorRetVA((ncor!=n) && (ncor!=n*n),hdf5_base,"Bad size for %s in %s (got %d expected %d or %d)",*err,__LINE__,,"cor",cur_lkl,ncor,n,n*n);
    isdiag = 0;
    if (ncor==n) {
      isdiag = 1;
    } 
  }
  
  *_al = al;
  *_bl = bl;
  *_nl = nl;
  *_cor = corr;
  *_isdiag = isdiag;
  
}

cmblkl* clik_ivg_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  tease *ing;
  cmblkl *cing;
  int i,isdiag;
  double *al, *bl, *nl,*corr;
  int nlcst;
  double zero;
  int n;
  int *Ml;
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  
  zero = 0;
  n=nell;
  Ml = ell;
  if (nbins!=0) {
    n=nbins;
    Ml = NULL;
  }
  
  al_bl_cor(group_id,cur_lkl, &n, &al, &bl, &nl, &corr, &isdiag,err);
  forwardError(*err,__LINE__,NULL);
  
  nlcst = 0;
  if (nl == NULL) {
    nl = &zero;
    nlcst = 1;    
  }
    
  ing = tease_init(n,Ml,al, bl, nl, nlcst, corr, isdiag, err);
  forwardError(*err,__LINE__,NULL);
  
  free(al);
  free(bl);
  if (nl!=&zero) {
    free(nl);
  }
  if (corr!=NULL) {
    free(corr);
  }
  
  cing = init_cmblkl(ing, &tease_log_pdf, 
                     &tease_free,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}


cmblkl* clik_gauss_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  gausslkl *ing;
  cmblkl *cing;
  int i,isdiag;
  double *al, *bl, *nl,*corr;
  int nlcst;
  double zero;
  int n;
  int *Ml;
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;  
  herr_t hstat;
  
  zero = 0;
  n=nell;
  Ml = ell;
  if (nbins!=0) {
    n=nbins;
    Ml = NULL;
  }
  
  al_bl_cor(group_id,cur_lkl, &n, &al, &bl, &nl, &corr, &isdiag,err);
  forwardError(*err,__LINE__,NULL);
  
  nlcst = 0;
  if (nl == NULL) {
    nl = &zero;
    nlcst = 1;    
  }
  
  ing = gausslkl_init(n,Ml,al, bl, nl, nlcst, corr, isdiag, err);
  forwardError(*err,__LINE__,NULL);
  
  free(al);
  free(bl);
  if (nl!=&zero) {
    free(nl);
  }
  if (corr!=NULL) {
    free(corr);
  }
  
  cing = init_cmblkl(ing, &gausslkl_log_pdf, 
                     &tease_free,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;  
}


cmblkl* clik_smica_init(hid_t group_id, char* cur_lkl, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  Smica *ing;
  cmblkl *cing;  
  double zero;
  int n;
  int *Ml;
  double *wq,*rq_hat,*rq_0;
  int m;
  SmicaComp **SCs;
  int ic;
  int mT,mP,nc;
  double *A_cmb;
  Smica *smic;
  int ncl,icl,nb;
  hsize_t ndum;
  H5T_class_t dum;
  size_t ddum;
  herr_t hstat;
  int cnt,xdim;
  char **xnames;
  parname *xnames_buf;
  
  zero = 0;
  testErrorRet(nbins==0,-101010,"no binning matrix. Argl",*err,__LINE__,NULL);
  ncl = 0;
  for(icl=0;icl<6;icl++) {
    if (has_cl[icl]==1) {
      ncl++;
    }
  }
  
  nb = nbins/ncl;
  testErrorRet(nbins!=nb*ncl,-101010,"bad binning matrix. Argl",*err,__LINE__,NULL);
  
  
  // try to read the bin weights
  wq = NULL;
  hstat = H5LTfind_dataset(group_id, "wq");
  if (hstat==1) {
    wq = hdf5_double_datarray(group_id, cur_lkl,"wq",&nb,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  
  // read the number of channels
  hstat = H5LTget_attribute_int( group_id, ".", "m_channel_T",  &mT);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read m_channel_T in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  hstat = H5LTget_attribute_int( group_id, ".", "m_channel_P",  &mP);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read m_channel_P in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  m = mT + mP;
  
  // read rq_hat
  int nrq;
  
  nrq = nb*m*m;
  rq_hat = hdf5_double_datarray(group_id, cur_lkl,"Rq_hat",&nrq,err);
  forwardError(*err,__LINE__,NULL);  
  
  // try to read rq_0
  rq_0 = NULL;
  hstat = H5LTfind_dataset(group_id, "Rq_0");
  if (hstat==1) {
    rq_0 = hdf5_double_datarray(group_id, cur_lkl,"Rq_0",&nrq,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  // how many components ?
  hstat = H5LTget_attribute_int( group_id, ".", "n_component",  &nc);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read n_component in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  
  SCs = malloc_err(sizeof(SmicaComp*) * nc,err);
  forwardError(*err,__LINE__,NULL);
  
  // now deal with the CMB component
  // read A_cmb
  A_cmb = hdf5_double_attarray(group_id,cur_lkl,"A_cmb",&m,err);
  forwardError(*err,__LINE__,NULL);    

  // init cmb comp
  SCs[0] = comp_CMB_init(nb, mT,mP, has_cl, A_cmb, err);
  forwardError(*err,__LINE__,NULL);    
  
  free(A_cmb);
  
  // deal with other components
  xdim = 0;
  for(ic=1;ic<nc;ic++) {
    char cur_cmp[256];
    char cur_cmp_tot[256];
    clik_smica_comp_init_func *smica_dl_init;
    void* dlhandle;
    char init_func_name[256];
    
#ifdef HAS_RTLD_DEFAULT 
    dlhandle = RTLD_DEFAULT;
#else
    dlhandle = NULL;
#endif
    parname comp_type;
    hid_t comp_id;
    
    SCs[ic] = NULL;
    
    sprintf(cur_cmp,"component_%d",ic);
    sprintf(cur_cmp_tot,"%s/component_%d",cur_lkl,ic);
    
    comp_id = H5Gopen(group_id, cur_cmp, H5P_DEFAULT);
    testErrorRetVA(comp_id<0,hdf5_base,"cannot read component %s in %s (got %d)",*err,__LINE__,NULL,cur_cmp,cur_lkl,hstat);
    
    // get type
    memset(comp_type,0,_pn_size*sizeof(char));
    hstat = H5LTget_attribute_string(comp_id, ".", "component_type",  comp_type);
    testErrorRetVA(hstat<0,hdf5_base,"cannot read component_type in %s/%s (got %d)",*err,__LINE__,NULL,cur_lkl,cur_cmp,hstat);
    
    
    sprintf(init_func_name,"clik_smica_comp_%s_init",comp_type);
    smica_dl_init = dlsym(dlhandle,init_func_name);
    testErrorRetVA(smica_dl_init==NULL,-1111,"Cannot initialize smica component type %s from %s dl error : %s",*err,__LINE__,NULL,comp_type,cur_lkl,dlerror()); 

    SCs[ic] = smica_dl_init(comp_id,cur_cmp_tot,nb,m, nell, ell, has_cl, unit, wl, bins,nbins,err);
    forwardError(*err,__LINE__,NULL);
    
    hstat = H5Gclose(comp_id);
    testErrorRetVA(hstat<0,hdf5_base,"cannot close %s in  %s (got %d)",*err,__LINE__,NULL,cur_cmp,cur_lkl,hstat);    
    
    xdim += SCs[ic]->ndim;
    
  }
  
  // deal with names and xdims
  if (xdim!=0) {
    xnames = malloc_err(sizeof(char*)*xdim,err);
    forwardError(*err,__LINE__,NULL);
  
    xnames_buf = malloc_err(sizeof(parname)*xdim,err);
    forwardError(*err,__LINE__,NULL);
    cnt = 0;
    for(ic=1;ic<nc;ic++) {
      int ix;
      for(ix=0;ix<SCs[ic]->ndim;ix++) {
        if (SCs[ic]->names!=NULL) {
          xnames[cnt]=&(SCs[ic]->names[ix]);
        } else {
          sprintf(xnames_buf[cnt],"SMICA_COMP_%d_%d",ic,ix);
          xnames[cnt]=&(xnames_buf[cnt]);
        }
        cnt++;
      }
    }
  }
  smic = Smica_init(nb, wq, m, rq_hat, rq_0, nc, SCs,err);
  forwardError(*err,__LINE__,NULL);
  
  free(rq_hat);
  
  if (rq_0!=NULL) {
    free(rq_0);
  }    
  
  if (wq!=NULL) {
    free(wq);
  }  
  
  cing = init_cmblkl(smic, &Smica_lkl, 
                     &free_Smica,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  if (xdim!=0) {
    cmblkl_set_names(cing, xnames,err);
    forwardError(*err,__LINE__,NULL);

    free(xnames);
    free(xnames_buf);
  }
  
  return cing;  
}

SmicaComp * clik_smica_comp_1d_init(hid_t comp_id, char* cur_lkl,int nb, int m,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  double *someA;
  int nA;
  SmicaComp *SC;
  herr_t hstat;
  
  // try to read A

  someA = NULL;
  hstat = H5LTfind_attribute(comp_id, "A");
  if (hstat==1) {
    nA = -1;
    someA = hdf5_double_attarray(comp_id,cur_lkl,"A",&nA,err);
  }
  testErrorRetVA(someA!=NULL && nA!=m,-100,"Not enough data in %s (expected %d got %d)",*err,__LINE__,NULL,cur_lkl,nA,m);
  SC = comp_1D_init(nb, m, someA, err);
  forwardError(*err,__LINE__,NULL);    
  if (someA!=NULL) {
    free(someA);
  }
  return SC;
}

SmicaComp * clik_smica_comp_nd_init(hid_t comp_id, char* cur_lkl,int nb, int m,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  double *someA;
  int nA;
  int nd;
  SmicaComp *SC;
  herr_t hstat;
  // try to read A

  hstat = H5LTget_attribute_int(comp_id, ".", "nd",  &nd);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read nd in component %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  someA = NULL;
  hstat = H5LTfind_attribute(comp_id, "A");
  if (hstat==1) {
    nA = -1;
    someA = hdf5_double_attarray(comp_id,cur_lkl,"A",&nA,err);
  }
  
  testErrorRetVA(someA!=NULL && nA!=m*nd,-100,"Not enough data in %s (expected %d got %d)",*err,__LINE__,NULL,cur_lkl,nA,m*nd);
  SC = comp_nD_init(nb, m, nd, someA, err);
  forwardError(*err,__LINE__,NULL);    
  if (someA!=NULL) {
    free(someA);
  }
  return SC;
}

SmicaComp * clik_smica_comp_diag_init(hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  herr_t hstat;
  // try to read A

  SC = comp_diag_init(nb, m, err);
  forwardError(*err,__LINE__,NULL);    
  return SC;
}

SmicaComp * clik_smica_comp_cst_init(hid_t comp_id, char* cur_lkl,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  herr_t hstat;
  double *rq_0;
  int tt;

  rq_0 =  hdf5_double_datarray(comp_id,cur_lkl,"Rq_0",&tt,err);
  forwardError(*err,__LINE__,NULL);    
  _DEBUGHERE_("%d %d %d %d",m,nb,nell,tt);
  testErrorRetVA(tt != m*m*nb,-22345,"%s:cst component does not have the correct number of data (expected %d got %d)",*err,__LINE__,NULL,cur_lkl,m*m*nb,tt)

  SC = comp_cst_init(nb, m, rq_0, err);
  forwardError(*err,__LINE__,NULL);    

  free(rq_0);
  return SC;
}

SmicaComp * clik_smica_gcal_log_init(hid_t comp_id, char* cur_lkl,int nb, int m,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int *ngcal;
  double *gcaltpl;
  SmicaComp *SC;
  herr_t hstat;
  int mm,im,tt,ig;
  char **xnames,*bnames;
  int binned;

  mm = m;
  ngcal = hdf5_int_attarray(comp_id,cur_lkl,"ngcal",&mm,err);
  forwardError(*err,__LINE__,NULL);    

  tt = 0;
  for(im=0;im<m;im++) {
    testErrorRetVA(ngcal[im]<0,hdf5_base,"%s: ngcal[%d] does make any sence (got %d)",*err,__LINE__,NULL,cur_lkl,im,ngcal[im]);
    tt += ngcal[im];
  }
  
  gcaltpl =  hdf5_double_datarray(comp_id,cur_lkl,"gcaltpl",&tt,err);
  forwardError(*err,__LINE__,NULL);    

  hstat = H5LTget_attribute_int(comp_id, ".", "binned",  &binned);
  testErrorRetVA(hstat<0,hdf5_base,"cannot read binned in component %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  if (binned!=0) {
    binned = 1;
  }

  SC = comp_gcal_log_init(nb,m, ngcal, gcaltpl,nell*binned,bins,err);
  forwardError(*err,__LINE__,NULL);    

  hstat = H5LTfind_attribute(comp_id, "names");
  if (hstat==1) {
    int dz;
    dz =-1;
    bnames = hdf5_char_attarray(comp_id,cur_lkl,"names",&dz, err);
    forwardError(*err,__LINE__,NULL); 
  } else {
    int ii;
    bnames = malloc_err(sizeof(char)*256*tt,err);
    forwardError(*err,__LINE__,NULL); 
    ii=0;
    for(im=0;im<m;im++) {
      for(ig=0;ig<ngcal[im];ig++) {
        sprintf(&(bnames[ii*256]),"gcal_%d_%d",im,ig);
        ii++;
      }
    }
  }

  xnames = malloc_err(sizeof(char*)*tt,err);
  for(im=0;im<tt;im++) {
    xnames[im] =&(bnames[im*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames); 
  free(bnames); 
  free(ngcal);
  free(gcaltpl);

  return SC;
}

SmicaComp * clik_smica_gcal_lin_init(hid_t comp_id, char* cur_lkl,int nb, int m,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int *ngcal;
  double *gcaltpl;
  SmicaComp *SC;
  herr_t hstat;
  int mm,im,tt;

  mm = m;
  ngcal = hdf5_int_attarray(comp_id,cur_lkl,"ngcal",&mm,err);
  forwardError(*err,__LINE__,NULL);    

  tt = 0;
  for(im=0;im<m;im++) {
    testErrorRetVA(ngcal[im]<0,hdf5_base,"%s: ngcal[%d] does make any sence (got %d)",*err,__LINE__,NULL,cur_lkl,im,ngcal[im]);
    tt += ngcal[im];
  }
  
  gcaltpl =  hdf5_double_datarray(comp_id,cur_lkl,"gcaltpl",&tt,err);
  forwardError(*err,__LINE__,NULL);    
  SC = comp_gcal_lin_init(nb,m, ngcal, gcaltpl,nell,bins,err);
  forwardError(*err,__LINE__,NULL);    

  free(ngcal);
  free(gcaltpl);

  return SC;
}
