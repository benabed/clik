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
     //grow list;
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
    //_DEBUGHERE_("%p",value);
    //_DEBUGHERE_("%p",value[0]);
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

  epl = malloc_err(sizeof(parametric),err);
  forwardError(*err,__LINE__,NULL);
  
  epl->pf = pflist_init(err);
  forwardError(*err,__LINE__,NULL);
  
  epl->default_settings = pflist_init(err);
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
  
  if (epl->nvar>0) {
    epl->sdRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq*epl->nvar,err);
    forwardError(*err,__LINE__,NULL);
  } else {
      epl->sdRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq*1,err);
      forwardError(*err,__LINE__,NULL);
  }
  
  epl->varkey = &(epl->pf->key[epl->ndef]);
  
  epl->dnofail=0;

  return epl;
}

void parametric_dnofail(parametric* egl, int vl) {
  egl->dnofail = vl;
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
  pflist_free(&(egl->pf));
  pflist_free(&(egl->default_settings));
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

void parametric_end_derivative_loop(parametric *egl,double* dRq, char* varkey, error **err) {
  testErrorRetVA(egl->dnofail!=1,-1234,"Cannot derive on parameter '%s'",*err,__LINE__,,varkey);
  memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*egl->nfreq*egl->nfreq);
}

// ASTRO PART //

double parametric_get_default(parametric* egl,char *key, error **err) {
  char *res;
  double vres;
  res = pflist_get_value(egl->default_settings,key,NULL,err);
  forwardError(*err,__LINE__,0);
  testErrorRetVA(res==NULL,-123432,"unknown default setting '%s'",*err,__LINE__,-1,key);
  //_DEBUGHERE_("'%s' -> '%s'",key,res);
  vres = atof(res);
  //_DEBUGHERE_("%g",vres);
  return vres;
}

void parametric_set_default(parametric* egl,char *key, double value,error **err) {
  char calue[200];
  char *nalue;

  sprintf(calue,"%30g",value);
  //_DEBUGHERE_("%p",egl->default_settings);
  //_DEBUGHERE_("%s",calue);
  nalue = calue;
  //_DEBUGHERE_("%p %p",&(nalue),calue);
  
  pflist_add_item(egl->default_settings,1, &key, (char**)&(nalue),err);
  forwardError(*err,__LINE__,);
  
  
}

double parametric_get_value(parametric *egl, char *key, error **err) {
  int ps;
  double res;

  ps=pflist_key_index(egl->pf,key,err);
  forwardError(*err,__LINE__,0);
  if (ps==-1) { // not in the default/var try in the default_settings
    res = parametric_get_default(egl,key,err);
    forwardError(*err,__LINE__,0);
    return res;
  }
  res = 0;
  res = pflist_get_double_value(egl->pf,key,&res,err);
  forwardError(*err,__LINE__,);
  return res;
}

// To get from intensity to \delta_T (CMB)

double dBdT(double nu, double nu0) {

  double x,x0,ex,ex0,res,res0;
  x0 = nu0/56.78;
  ex0 = exp(x0);
  x = nu/56.78;
  ex = exp(x);
  res0 = pow(x0,4.0)*ex0/((ex0-1.)*(ex0-1.));
  res = pow(x,4.0)*ex/((ex-1.)*(ex-1.));
  return(res/res0);

}

// Simple power low in ell, constant emissivity

parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_compute;
  egl->eg_free = NULL;
  
  // set default settings (those will be overide by defkey/value and varkey/value)
  parametric_set_default(egl,"l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"index",0,err);
  forwardError(*err,__LINE__,NULL);
   
  parametric_set_default(egl,"A",1,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

void powerlaw_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,A,v;
  
  egl = exg;
  l_pivot = parametric_get_value(egl,"l_pivot",err);
  forwardError(*err,__LINE__,);
  
  index = parametric_get_value(egl,"index",err);
  forwardError(*err,__LINE__,);

  A = parametric_get_value(egl,"A",err);  
  forwardError(*err,__LINE__,);

  
  nfreq = egl->nfreq;
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    //_DEBUGHERE_("%g %g %g %g",A,ell,l_pivot,index);
    v = A*pow((double) ell/l_pivot,(double) index);
    //_DEBUGHERE_("%g ",v);
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
          v = A * log((double)ell/l_pivot) * pow((double) ell/l_pivot,(double) index);
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
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }

  return;
}

// Simple power law in ell, arbitrary emissivity (including arbitrary cross-correlations, i.e. NOT rank=1 a priori)
// This could be used e.g for CIB Poisson

parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_free_emissivity_compute;
  egl->eg_free = &powerlaw_free_emissivity_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"index",0,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"A_%d_%d",egl->freqlist[m1],egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

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
  l_pivot = parametric_get_value(egl,"l_pivot",err);
  forwardError(*err,__LINE__,);
  
  index = parametric_get_value(egl,"index",err);
  forwardError(*err,__LINE__,);

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"A_%d_%d",egl->freqlist[m1],egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
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
          v = log(ell/l_pivot) * pow(ell/l_pivot,index);
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
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }
  return;
}

// Millea et al. radio galaxies Poisson contribution. Based on a single population with given mean spectral index, as well as dispersion

parametric *radiogal_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &radiogal_compute;
  egl->eg_free = &radiogal_free;
  egl->payload = malloc_err(sizeof(double)*(2*egl->nfreq*egl->nfreq + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"radiogal_norm",78.5,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"radiogal_alpha",-0.36,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"radiogal_sigma",0.64,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void radiogal_compute(void* exg, double *Rq, double* dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm_rg, alpha_rg, sigma_rg;
  double d3000, nu0;
  double x, lnx, ln2x;
  double *A,*B,*vec;

  egl = exg;
  A = egl->payload;
  d3000 = 3000.*3001./2./M_PI;
  nu0 = 143;

  norm_rg = parametric_get_value(egl,"radiogal_norm",err);
  forwardError(*err,__LINE__,);

  alpha_rg = parametric_get_value(egl,"radiogal_alpha",err);
  forwardError(*err,__LINE__,);

  sigma_rg = parametric_get_value(egl,"radiogal_sigma",err);
  forwardError(*err,__LINE__,);

  nfreq = egl->nfreq;
  B = A + nfreq*nfreq;
  vec = B + nfreq*nfreq;
  for (m1=0;m1<nfreq;m1++) {
    vec[m1] = dBdT((double)egl->freqlist[m1],nu0);
    for(m2=m1;m2<nfreq;m2++) {
      x = (double)egl->freqlist[m1]*(double)egl->freqlist[m2]/(nu0*nu0);
      lnx = log(x);
      ln2x = lnx*lnx;
      A[m1*nfreq+m2] = lnx;
      A[m2*nfreq+m1] = lnx;
      B[m1*nfreq+m2] = ln2x;
      B[m2*nfreq+m1] = ln2x;
    }
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    mell=(ell-egl->lmin)*nfreq*nfreq;
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
        Rq[mell + m1*nfreq+m2] = norm_rg/d3000 * 
          exp(alpha_rg*A[m1*nfreq+m2] +
              sigma_rg*sigma_rg/2.0 * B[m1*nfreq+m2])/(vec[m1]*vec[m2]);
        Rq[mell + m2*nfreq+m1] = Rq[mell + m1*nfreq+m2];
      }
    }
  }

  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      if (strcmp(egl->varkey[iv],"radiogal_norm")==0) {
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = Rq[mell+m1*nfreq+m2]/norm_rg;
              dRq[mv+mell + m2*nfreq+m1] = Rq[mell+m2*nfreq+m1]/norm_rg;
            }
          }
        }
        continue;
      }
      if (strcmp(egl->varkey[iv],"radiogal_alpha")==0) {
        // dR/dalpha_rg = log(nu1*nu2/nu0^2) * R
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = A[m1*nfreq+m2] * Rq[mell+m1*nfreq+m2];
              dRq[mv+mell + m2*nfreq+m1] = A[m2*nfreq+m1] * Rq[mell+m2*nfreq+m1];
            }
          }
        }
        continue;
      }
      if (strcmp(egl->varkey[iv],"radiogal_sigma")==0) {
        // dR/dsigma_rg = log(nu1*nu2/nu0^2)^2 * R
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              dRq[mv+mell + m1*nfreq+m2] = sigma_rg * B[m1*nfreq+m2] * Rq[mell+m1*nfreq+m2];
              dRq[mv+mell + m2*nfreq+m1] = sigma_rg * B[m2*nfreq+m1] * Rq[mell+m2*nfreq+m1];
            }
          }
        }
        continue;
      }
      // error return
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }
  return;
}

void radiogal_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

// Galactic model

// Dust emissivity normalized to 1 at nu0
double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res, res0;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  res0 = pow(nu0,3.0+beta_dust)/(ex0 - 1.); // Overall normalization will go away
  x = nu * h_over_kT;
  ex = exp(x);
  res = pow(nu,3.0+beta_dust)/(ex - 1.);

  // Return dust emissivity normalized to one at nu0, in dT (CMB)
  return (res/res0/dBdT(nu,nu0));
}

// Derivative of dust_spectrum with respect to T_dust
double d_dust_spectrum_d_T_dust(double nu, double T_dust, double beta_dust, double nu0) {

  double h_over_kT, x, x0, ex, ex0, res;

  h_over_kT = 1./(T_dust * 20.836); // Assumes frequencies in GHz
  x0 = nu0 * h_over_kT;
  ex0 = exp(x0);
  x = nu * h_over_kT;
  ex = exp(x);
  res = 1./T_dust * pow(nu/nu0,3.0+beta_dust) * (ex0-1.) * x/((ex - 1.)*(ex - 1.));
  res /= dBdT(nu,nu0);

  return(res);
  
} 

// Derivative of dust_spectrum with respect to beta_dust
double d_dust_spectrum_d_beta_dust(double nu, double T_dust, double beta_dust, double nu0) {
  
  return (log(nu/nu0)*dust_spectrum(nu,T_dust,beta_dust,nu0));

}

// Non-thermal emissivity normalized to 1 at nu0 (e.g. synchrotron, free-free, etc.)
// The spectral index is called alpha here, to distinguish from the gray body case
// Expressed in dT (CMB)
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0) {

  return (pow(nu/nu0,alpha_non_thermal)/dBdT(nu,nu0));
}

// Derivative of non_thermal_spectrum with respect to alpha
double d_non_thermal_spectrum_d_alpha_non_thermal(double nu, double alpha_non_thermal, double no) {

  return (log(nu/nu0)*non_thermal_spectrum(nu,alpha_non_thermal,nu0));
}

parametric *galactic_component_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  pfchar type;
  char* pt;
  int isdust;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &galactic_component_compute;
  egl->eg_free = &galactic_component_free;
  egl->payload = malloc_err(sizeof(double)*2*egl->nfreq*(egl->nfreq+1),err);
  forwardError(*err,__LINE__,NULL);
  
  // uK^2 at l=500, nu=143 GHz;
  parametric_set_default(egl,"gal_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal_index",0,err);
  forwardError(*err,__LINE__,NULL);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,);

  if (strcmp(type,"dust")==0) {
    parametric_set_default(egl,"gal_beta_dust",1.8,err);
    forwardError(*err,__LINE__,NULL);

    parametric_set_default(egl,"gal_T_dust",1.8,err);
    forwardError(*err,__LINE__,NULL);
  } else if (strcmp(type,"non_thermal")==0) {
    //Intensity, = -3.0 in RJ
    parametric_set_default(egl,"gal_alpha_non_thermal",-1,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,,type);
    // return ?
  }

  return egl;
}

void galactic_component_compute(void* exg, double *Rq, double *dRq, error **err) {
  parametric *egl;
  int ell,m1,m2,mell,nfreq,iv,mv;
  double norm,l_pivot,index,alpha_non_thermal,beta_dust,T_dust;
  double v,lA,nu0;
  double *A, *a, *B, *b;
  pfchar type;
  char* pt;
  int isdust;

  nu0 = 143;
  egl = exg;

  // dust or non_thermal

  norm = parametric_get_value(egl,"gal_norm",err);
  forwardError(*err,__LINE__,);

  l_pivot = parametric_get_value(egl,"gal_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  index = parametric_get_value(egl,"gal_index",err);
  forwardError(*err,__LINE__,);

  sprintf(type,"dust");
  isdust = 1;
  pt = (char*)type;
  pt = pflist_get_value(egl->pf,"gal_type",pt,err);
  forwardError(*err,__LINE__,);

  if (strcmp(type,"dust")==0) {
    
    beta_dust = parametric_get_value(egl,"gal_beta_dust",err);
    forwardError(*err,__LINE__,);

    T_dust = parametric_get_value(egl,"gal_T_dust",err);
    forwardError(*err,__LINE__,);

  } else if (strcmp(type,"non_thermal")==0) {

    isdust = 0;
    alpha_non_thermal = parametric_get_value(egl,"gal_alpha_non_thermal",err);
    forwardError(*err,__LINE__,);

  } else {
    testErrorRetVA(1==1,-1234,"Unknown Galactic component type '%s'",*err,__LINE__,,type);
  }
    

  A = egl->payload;            // Rank-one, tensorial emissivity
  nfreq = egl->nfreq;
  a = &(A[nfreq*nfreq]);       // Emissivity vector

  B = &(A[nfreq*(nfreq+1)]);   // Used for derivatives
  b = &(A[nfreq*(2*nfreq+1)]); // Idem

  for (m1=0;m1<nfreq;m1++) {
    if (isdust) {
      a[m1] = dust_spectrum((double)egl->freqlist[m1],T_dust,beta_dust,nu0);
    } else {
      a[m1] = non_thermal_spectrum((double)egl->freqlist[m1],alpha_non_thermal,nu0);
    }
  }

  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      v = a[m1]*a[m2];
      A[m1*nfreq+m2] = v;
      A[m2*nfreq+m1] = v;
    }
  }

  // R_q

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index) * norm;
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[mell + m1*nfreq + m2] = lA*v;
        Rq[mell + m2*nfreq + m1] = lA*v;
      }  
    }
  }

  // dR_q
  if (dRq!=NULL) {
    for (iv=0;iv<egl->nvar;iv++) {
      mv = iv*(egl->lmax-egl->lmin+1)*nfreq*nfreq;

      if (strcmp(egl->varkey[iv],"gal_norm")==0) {
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = pow(ell/l_pivot,index);
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              lA = A[m1*nfreq+m2];
              dRq[mv+mell + m1*nfreq+m2] = lA*v;
              dRq[mv+mell + m2*nfreq+m1] = lA*v;
            }
          }
        }
        continue;
      }

      if (strcmp(egl->varkey[iv],"gal_index")==0) {
        // dR/dindex
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = log(ell/l_pivot)*pow(ell/l_pivot,index) * norm;
          mell = (ell-egl->lmin)*nfreq*nfreq;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              lA = A[m1*nfreq+m2];
              dRq[mv+mell + m1*nfreq+m2] = lA*v;
              dRq[mv+mell + m2*nfreq+m1] = lA*v;
            }
          }
        }
        continue;
      }

      if (strcmp(egl->varkey[iv],"gal_beta_dust")==0 && isdust) {
        // dR/dbeta_dust
        // Get vector emissivity derivative
        for (m1=0;m1<nfreq;m1++) {
          b[m1] = d_dust_spectrum_d_beta_dust((double)egl->freqlist[m1],T_dust,beta_dust,nu0);
          for (m2=m1;m2<nfreq;m2++) {
            B[m1*nfreq+m2] = b[m1]*b[m2];
            B[m2*nfreq+m1] = b[m1]*b[m2];
          }
        }
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          v = power(ell/l_pivot,index) * norm;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              lA = B[m1*nfreq+m2];
              dRq[mv+mell + m1*nfreq+m2] = lA*v;
              dRq[mv+mell + m2*nfreq+m1] = lA*v;
            }
          }
        }
        continue;
      }

      if (strcmp(egl->varkey[iv],"gal_T_dust")==0 && isdust) {
        // dR/dT_dust
        // Get vector emissivity derivative
        for (m1=0;m1<nfreq;m1++) {
          b[m1] = d_dust_spectrum_d_T_dust((double)egl->freqlist[m1],T_dust,beta_dust,nu0);
          for (m2=m1;m2<nfreq;m2++) {
            B[m1*nfreq+m2] = b[m1]*b[m2];
            B[m2*nfreq+m1] = b[m1]*b[m2];
          }
        }
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          v = power(ell/l_pivot,index) * norm;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              lA = B[m1*nfreq+m2];
              dRq[mv+mell + m1*nfreq+m2] = lA*v;
              dRq[mv+mell + m2*nfreq+m1] = lA*v;
            }
          }
        }
        continue;
      }

      if (strcmp(egl->varkey[iv],"gal_alpha_non_thermal")==0 && (isdust==0)) {
        // dR/dalpha_non_thermal
        // Get vector emissivity derivative
        for (m1=0;m1<nfreq;m1++) {
          b[m1] = d_non_thermal_spectrum_d_alpha_non_thermal((double)egl->freqlist[m1],alpha_non_thermal,nu0);
          for (m2=m1;m2<nfreq;m2++) {
            B[m1*nfreq+m2] = b[m1]*b[m2];
            B[m2*nfreq+m1] = b[m1]*b[m2];
          }
        }
        for (ell=egl->lmin;ell<=egl->lmax;ell++) {
          mell = (ell-egl->lmin)*nfreq*nfreq;
          v = power(ell/l_pivot,index) * norm;
          for (m1=0;m1<nfreq;m1++) {
            for (m2=m1;m2<nfreq;m2++) {
              lA = B[m1*nfreq+m2];
              dRq[mv+mell + m1*nfreq+m2] = lA*v;
              dRq[mv+mell + m2*nfreq+m1] = lA*v;
            }
          }
        }
        continue;
      }

      // error return
      parametric_end_derivative_loop(egl,&(dRq[mv]),egl->varkey[iv],err);
      forwardError(*err,__LINE__,);
    }
  }
  
}

void galactic_component_free(void **pp) {
  double *p;

  p = *pp;
  free(p);
  *pp=NULL;
}

