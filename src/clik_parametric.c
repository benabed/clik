#include "clik_parametric.h"


pflist* pflist_init(error **err) {
 pflist *pf;
 
 pf = malloc_err(sizeof(pflist),err);
 forwardError(*err,__LINE__,NULL);
 
 pf->nmax = 20;
 pf->nkey = 0;
 pf->key = malloc_err(sizeof(pfchar)*pf->nmax,err);
 forwardError(*err,__LINE__,NULL);
 pf->value = malloc_err(sizeof(pfchar)*pf->nmax,err);
 forwardError(*err,__LINE__,NULL);
    
return pf;
}

void pflist_free(void **ppf) {
  pflist *pf;

  pf = *ppf;
  if (pf->nmax !=0) {
    free(pf->key);
    free(pf->value);
  }
  free(pf);
  *ppf = NULL;
}

void pflist_add_item(pflist* pf,int nit, char** key, char **value,error **err) {
  int i;
  int cur;

  if (pf->nkey+nit>= pf->nmax) {
    // grow list;
    //_DEBUGHERE_("-> %p %d %d %d",pf,pf->nmax, pf->nkey, nit);

    int nnmax =  (pf->nkey+nit)*2;
    pf->key=resize_err(pf->key, sizeof(pfchar)*pf->nmax, sizeof(pfchar)*nnmax, 1, err);
    forwardError(*err,__LINE__,);
    pf->value=resize_err(pf->value, sizeof(pfchar)*pf->nmax, sizeof(pfchar)*nnmax, 1, err);
    forwardError(*err,__LINE__,);
    pf->nmax =  nnmax;
  }
  //_DEBUGHERE_("-> %p %d %d %d",pf,pf->nmax, pf->nkey, nit);
  cur = pf->nkey;
  for(i=0;i<nit;i++) {
    int idx;
    idx = pflist_key_index(pf,key[i],err);
    forwardError(*err,__LINE__,);
    //_DEBUGHERE_("cur %d idx %d",cur,idx)
    if (idx==-1) {
      idx = cur;
      cur++;
    }
    //_DEBUGHERE_("cur %d idx %d",cur,idx)
    strcpy(pf->key[idx],key[i]);
    if (value[i] !=NULL) {
      strcpy(pf->value[idx],value[i]);  
    } else {
      pf->value[idx][0] = '\0';
    }
    //_DEBUGHERE_("added %d '%s' : '%s'",idx,pf->key[idx],pf->value[idx]);
  }
  pf->nkey = cur;  
  //_DEBUGHERE_("-> %p %d %d %d",pf,pf->nmax, pf->nkey, nit);

}

int pflist_key_index(pflist *pf, char *key, error **err) {
  int i;

  for(i=0;i<pf->nkey;i++) {
    if (strcmp(key,pf->key[i])==0) {
      return i;
    }
  }
  return -1;
}

char* pflist_get_value(pflist* pf, char* key, char* safeguard,error **err) {
  int ps;
  
  ps=pflist_key_index(pf,key,err);
  forwardError(*err,__LINE__,NULL);
  if (ps==-1) {
    testErrorRetVA(safeguard==NULL,-1234,"key '%s' absent",*err,__LINE__,NULL,key);
    return safeguard;
  }
  return pf->value[ps];
}

long pflist_get_int_value(pflist *pf, char *key,long* safeguard, error **err) {
  char *res;

  res = pflist_get_value(pf,key,(char*) safeguard,err);
  forwardError(*err,__LINE__,0);
  
  if (res==(char*)safeguard) {
    return *safeguard;
  }

  return atol(res);
}

double pflist_get_double_value(pflist *pf, char *key,double *safeguard, error **err) {
    char *res;

  res = pflist_get_value(pf,key,(char*)safeguard,err);
  forwardError(*err,__LINE__,0);
  
  if (res==(char*)safeguard) {
    return *safeguard;
  }

  return atof(res);
}

void pflist_remove_item(pflist* pf, int index,error **err) {
  int i;

  for(i=index;i<pf->nkey-1;i++) {
    strcpy(pf->key[i],pf->key[i+1]);
    strcpy(pf->value[i],pf->value[i+1]);
  }
  pf->nkey--;
}

void get_freq(int ndet, int *detlist, int* pnfreq, int** pfreqlist, int** pdet2freq, error **err) {
  int i,j;
  int *freqlist,*det2freq;
  int nfreq;
  
  freqlist = malloc_err(sizeof(int)*ndet,err);
  forwardError(*err,__LINE__,);
  
  det2freq = malloc_err(sizeof(int)*ndet, err);
  forwardError(*err,__LINE__,);
  
  nfreq = 0;
  for(i=0;i<ndet;i++) {
    det2freq[i]=-1;
    for(j=0;j<nfreq;j++) {
      if (freqlist[j]==detlist[i]) {
        det2freq[i] = j;
        break;
      }
    }
    if (det2freq[i]==-1) {
      det2freq[i]=nfreq;
      freqlist[nfreq]=detlist[i];
      nfreq++;
    }
  }
  *pnfreq = nfreq;
  *pfreqlist = freqlist;
  *pdet2freq = det2freq;
}


parametric *parametric_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *epl;
  int i;
  char *nop;

  //_DEBUGHERE_("","");
  epl = malloc_err(sizeof(parametric),err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("-> %p",epl);

  //_DEBUGHERE_("","");
  epl->pf = pflist_init(err);
  forwardError(*err,__LINE__,NULL);
  
  pflist_add_item(epl->pf,ndef,defkey,defvalue,err);
  forwardError(*err,__LINE__,NULL);
  
  nop = NULL;
  for(i=0;i<nvar;i++) {
    pflist_add_item(epl->pf,1,&(varkey[i]),&nop,err);
    forwardError(*err,__LINE__,NULL);
  }

  epl->nvar = nvar;
  epl->ndef = ndef;

  epl->ndet = ndet;
  get_freq(ndet, detlist, &(epl->nfreq), &(epl->freqlist), &(epl->det2freq),err);
  forwardError(*err,__LINE__,NULL);
  
  epl->payload = NULL;
  epl->eg_compute = NULL;
  epl->eg_free = NULL;

  epl->lmin = lmin;
  epl->lmax = lmax;

  epl->sRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  
  epl->sdRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq*epl->nvar,err);
  forwardError(*err,__LINE__,NULL);
  
  epl->varkey = &(epl->pf->key[epl->ndef]);
  
  return epl;
}

void parametric_free(void** pegl) {
  parametric *egl;

  egl = *pegl;
  if(egl->eg_free!=NULL) {
    egl->eg_free(&(egl->payload));
  }
    
  free(egl->det2freq);
    free(egl->freqlist);
    free(egl->sdRq);
    free(egl->sRq);
    free(egl);
    *pegl = NULL;
}

void parametric_compute(parametric *egl, double *pars, double* Rq, double *dRq, error **err) {
  int idet,jdet;
  int ell,dvar;
  int midet, mjdet;
  int m2,f2,mell,fell,ivar,mdvar,fdvar; 
  int dvaro;
  double *sRq, *sdRq;
  double v;

  
  testErrorRet(egl->eg_compute==NULL,-1234,"badly initialized",*err,__LINE__,);
  
  for(ivar=0;ivar<egl->nvar;ivar++) {
    sprintf(egl->pf->value[egl->ndef+ivar],"%40g",pars[ivar]);
  }
  
  if (egl->ndet!=egl->nfreq) {  
    sRq = NULL;
    if (Rq!=NULL) {
      sRq = egl->sRq;
    }
    sdRq = NULL;
    if (dRq!=NULL) {
      sdRq = egl->sdRq;
    }
    egl->eg_compute(egl,sRq, sdRq, err);
    forwardError(*err,__LINE__,);
  
    m2 = egl->ndet;
    m2 = m2*m2;
    f2 = egl->nfreq;
    f2 = f2*f2;
    mdvar = (egl->lmax+1-egl->lmin) * m2;
    fdvar = (egl->lmax+1-egl->lmin) * f2;

    if (sRq!=NULL) {
      for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
        mell = ell*m2;
        fell = ell*f2;
        for(idet=0;idet<egl->ndet;idet++) {
          midet = egl->det2freq[idet];
          for(jdet=idet;jdet<egl->ndet;jdet++) {
            mjdet = egl->det2freq[jdet];
            //_DEBUGHERE_("%d->%d %d->%d %d %d",idet,midet,jdet,mjdet,egl->ndet,egl->nfreq);
            v = sRq[fell + midet*egl->nfreq + mjdet];
            Rq[mell + idet*egl->ndet + jdet] = v;
            Rq[mell + jdet*egl->ndet + idet] = v;
          }
        }
      }      
    }
    if (sdRq!=NULL) {
      for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
        mell = ell*m2;
        fell = ell*f2; 
        for(idet=0;idet<egl->ndet;idet++) {
          midet = egl->det2freq[idet];
          for(jdet=idet;jdet<egl->ndet;jdet++) {
            mjdet = egl->det2freq[jdet];
            for(dvar=0;dvar<egl->nvar;dvar++) {
              v = sdRq[dvar*fdvar + fell + midet*egl->nfreq + mjdet];
              dRq[dvar*mdvar + mell + idet*egl->ndet + jdet] = v;
              dRq[dvar*mdvar + mell + jdet*egl->ndet + idet] = v;
            }
          }
        }
      }      
    }
  } else {
    egl->eg_compute(egl, Rq, dRq, err);
    forwardError(*err,__LINE__,);
  }
}


parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_compute;
  egl->eg_free = NULL;
  
  return egl;
}

void powerlaw_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,A,v;
  
  egl = exg;
  l_pivot = 500;
  l_pivot = pflist_get_double_value(egl->pf,"l_pivot",&l_pivot,err);
  forwardError(*err,__LINE__,);
  
  index = 0;
  index = pflist_get_double_value(egl->pf,"index",&index,err);
  forwardError(*err,__LINE__,);

  A = 1;
  A = pflist_get_double_value(egl->pf,"A",&A,err);
  forwardError(*err,__LINE__,);

  
  nfreq = egl->nfreq;
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = A*pow((double) ell/l_pivot,(double) index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[mell + m1*nfreq + m2] = v;
        Rq[mell + m2*nfreq + m1] = v;
      }  
    }
  }
  
  if (dRq!=NULL) {
    for(iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax+1-egl->lmin)*nfreq*nfreq;
      
      if (strcmp(egl->varkey[iv],"index")==0) {
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = index*A*pow((double) ell/l_pivot,(double) index-1);
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for(m1=0;m1<nfreq;m1++) {
            for(m2=m1;m2<nfreq;m2++) {
              //_DEBUGHERE_("%d %d %g",mv+mell + m1*nfreq + m2,mv+mell + m1*nfreq + m2,v)
              dRq[mv+mell + m1*nfreq + m2] = v;
              dRq[mv+mell + m2*nfreq + m1] = v;
            }  
          }
        }        
        continue;
      }
      
      if (strcmp(egl->varkey[iv],"A")==0) {
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = pow((double) ell/l_pivot,(double) index);
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for(m1=0;m1<nfreq;m1++) {
            for(m2=m1;m2<nfreq;m2++) {
              //_DEBUGHERE_("%d %d %g",mv+mell + m1*nfreq + m2,mv+mell + m1*nfreq + m2,v)
              dRq[mv+mell + m1*nfreq + m2] = v;
              dRq[mv+mell + m2*nfreq + m1] = v;
            }  
          }
        }        
        continue;
      }

      // error return
      testErrorRetVA(1==1,-1234,"Cannot derive on parameter '%s'",*err,__LINE__,,egl->varkey[iv]);
    }
  }

  return;
}

parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_free_emissivity_compute;
  egl->eg_free = &powerlaw_free_emissivity_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void powerlaw_free_emissivity_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

void powerlaw_free_emissivity_compute(void* exg, double *Rq, double *dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;

  egl = exg;
  l_pivot = 500;
  l_pivot = pflist_get_double_value(egl->pf,"l_pivot",&l_pivot,err);
  forwardError(*err,__LINE__,);
  
  index = 0;
  index = pflist_get_double_value(egl->pf,"index",&index,err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"A_%d_%d",egl->freqlist[m1],egl->freqlist[m2]);
      v = 1;
      v = pflist_get_double_value(egl->pf,name,&v,err);
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }


  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[mell + m1*nfreq + m2] = lA*v;
        Rq[mell + m2*nfreq + m1] = lA*v;
      }  
    }
  }

  if (dRq!=NULL) {
    for(iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax+1-egl->lmin)*nfreq*nfreq;
      
      if (strcmp(egl->varkey[iv],"index")==0) {
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = index*pow(ell/l_pivot,index-1);
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for(m1=0;m1<nfreq;m1++) {
            for(m2=m1;m2<nfreq;m2++) {
              lA = A[m1*nfreq+m2];
              //_DEBUGHERE_("%d %d %g %g %d %d",m1,m2,lA,v,mv+mell + m1*nfreq + m2,mv+mell + m2*nfreq + m1);
              dRq[mv+mell + m1*nfreq + m2] = lA*v;
              dRq[mv+mell + m2*nfreq + m1] = lA*v;
            }  
          }
        }        
        continue;
      }
      
      stop = 0;
      for(m1=0;m1<nfreq;m1++) {
        for(m2=m1;m2<nfreq;m2++) {
          sprintf(name,"A_%d_%d",egl->freqlist[m1],egl->freqlist[m2]);
  
          if (strcmp(egl->varkey[iv],name)==0) {
            stop=1;
            memset(&(dRq[mv]),0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
            for(ell=egl->lmin;ell<=egl->lmax;ell++) {
              v = pow(ell/l_pivot,index);
              mell = (ell-egl->lmin)*nfreq*nfreq;
              //_DEBUGHERE_("%d %d %g %g %d %d",m1,m2,lA,v,mv+mell + m1*nfreq + m2,mv+mell + m2*nfreq + m1);
              dRq[mv+mell + m1*nfreq + m2] = v;
              dRq[mv+mell + m2*nfreq + m1] = v;
            }
            break;
          }
        }
        if (stop==1) {
          break;
        }
      }
      if (stop==1) {
        continue;
      }
      
      // error return
      testErrorRetVA(1==1,-1234,"Cannot derive on parameter '%s'",*err,__LINE__,,egl->varkey[iv]);
    }
  }
  return;
}

