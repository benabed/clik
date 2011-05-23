/*
 *  distribution.c
 * simplified version of the pmclib full version
 */

#include "pmc.h"


/******************************************************************************/
/************** Distribution interface (should be moved somewhere else) *******/
/******************************************************************************/
distribution* init_simple_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                error **err) {
  distribution *dist;
  
  dist = init_distribution(ndim,data,log_pdf,freef,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  return dist;
}                                
                                  

distribution* init_distribution_full(int ndim,
                                     void* data, 
                                     posterior_log_pdf_func* log_pdf,
                                     posterior_log_free* freef,
                                     simulate_func *simulate,
                                     int n_ded,
                                     retrieve_ded_func* retrieve,
                                     error **err) {
  distribution *dist;
  
  testErrorRet(n_ded!=0 && retrieve==NULL,dist_undef,
               "Target invalid, expect deduced parameters, but no function provided...",
               *err,__LINE__,NULL);
  
  dist = malloc_err(sizeof(distribution),err);
  forwardError(*err,__LINE__,NULL);
  
  dist->ndim = ndim;
  dist->n_ded = n_ded;
  dist->data = data;
  dist->log_pdf = log_pdf;
  dist->free = freef;
  dist->simulate = simulate;
  dist->retrieve = retrieve;
  dist->broadcast_mpi = NULL;
  dist->ndef = 0;
  dist->def = NULL;
  dist->pars = NULL;
  dist->dlhandle = NULL;
  dist->name = NULL;
  
  dist->f_der = NULL;
  dist->d_der = NULL;
  return dist;
}

void distribution_set_broadcast(distribution* dist, mpi_exchange_func* broadcast, error **err) {
  testErrorRet(dist==NULL,-1,"Bad distribution",*err,__LINE__,);
  dist->broadcast_mpi = broadcast;
  return;
}

distribution* init_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                simulate_func *simulate,
                                error **err) {
  distribution  *dist;
  dist = init_distribution_full(ndim, data, log_pdf, freef, simulate, 0, NULL, err);
  forwardError(*err,__LINE__,NULL);
  return dist;
}

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err) {
  int i;
  
  testErrorRetVA(ndef >= dist->ndim + dist->ndef,dist_undef,"Too many defaults ! (expected at much %d, got %d)",*err,__LINE__,,dist->ndim+dist->ndef,ndef);

  if (dist->ndef!=0) {
    //discard old defaults
    free(dist->def);
    free(dist->pars);
    dist->ndim += dist->ndef;
    dist->ndef = 0;
  }
  
  dist->def = malloc_err(sizeof(int)*dist->ndim,err);
  forwardError(*err,__LINE__,);
  
  memset(dist->def,0,sizeof(int)*dist->ndim);
  
  dist->pars = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<ndef;i++) {
    testErrorRetVA(idef[i]>=dist->ndim,dist_undef,"Too many defaults ! (expected at much %d, got %d)",*err,__LINE__,,dist->ndim,idef[i]);
    testErrorRetVA(dist->def[idef[i]]==1,dist_undef,"par %d already set",*err,__LINE__,,idef[i]);
    
    dist->def[idef[i]] = 1;
    dist->pars[idef[i]] = vdef[i]; 
  }
  
  dist->ndef = ndef;
  dist->ndim -= ndef;
  
}

void distribution_set_default_name(distribution *dist, int ndef, char** namedef, double* vdef,error **err) {
  int* idef;
  
  idef = distribution_get_names(dist,ndef,namedef,0,err);
  forwardError(*err,__LINE__,);
  distribution_set_default(dist,ndef,idef,vdef,err);
  forwardError(*err,__LINE__,);
  free(idef);
}
int distribution_get_name(distribution *dist,char* name,error **err) {
  int i,j;
  
  testErrorRet(dist->name==NULL,dist_undef,"The distribution has no parameter name defined",*err,__LINE__,-1);
  j=0;
  for (i=0;i<dist->ndim + dist->ndef;i++) {
    if (dist->ndef>0 && dist->def[i]==1) {
      continue;
    } 
    if (strcmp(name,dist->name[i])==0) {
      return j;
    }
    j++;
  }

  j=-1;
  for (i=0;i<dist->n_ded;i++) {
    if (strcmp(name,dist->name[i])==0) {
      return j;
    }
    j--;
  }

  *err = addErrorVA(dist_undef, "Cannot find parameter %s", *err,__LINE__,name);
  return -1;
}

int* distribution_get_names(distribution *dist, int nnames, char** name, int includeded, error **err) {
  int* idef;
  int i;
  
  idef = malloc_err(sizeof(int)*nnames,err);
  forwardError(*err,__LINE__, NULL);
  
  for(i=0;i<nnames;i++) {
    idef[i] = distribution_get_name(dist,name[i],err);
    forwardError(*err,__LINE__, NULL);
    testErrorRetVA(idef[i]<0 && includeded==0,dist_undef,"Asking for deduced parameter (%s)",*err,__LINE__,NULL,name[i]);
  }
  
  return idef;
}


void distribution_set_names(distribution *dist,char** name, error **err) {
  int i;
  
  testErrorRet(dist->ndef!=0,dist_undef,"cannot set name after having set default parameters",*err,__LINE__,);
  if (dist->name!=NULL) {
    free(dist->name);
  }
  
  dist->name = malloc_err(sizeof(_char_name)*(dist->ndim+dist->n_ded),err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<dist->ndim+dist->n_ded;i++) {    
    sprintf(dist->name[i],"%s",name[i]);
  }
  
  return;
}

double * distribution_fill_pars(distribution *dist, double* pars, error **err) {
  const double *_pars;
  _pars = pars;
  if (dist->ndef!=0) {
    // i need to reset the parameters !
    int ir,ik;
    ir=0;
    for(ik=0;ik<dist->ndim+dist->ndef;ik++) {
      if (dist->def[ik]==0) {
        dist->pars[ik]=pars[ir];
        ir++;
      }
    }
    _pars = dist->pars;
  }
  return _pars;
}

double distribution_lkl(void* pdist, const double* pars, error **err) {
  distribution *dist;
  double res;
  const double *_pars;
  
  dist = pdist;
  testErrorRet(dist->log_pdf==NULL,dist_undef,"undefined log pdf function for distribution",*err,__LINE__,0);

  _pars = distribution_fill_pars(pdist, pars, err);
  forwardError(*err,__LINE__,0);
  
  res = dist->log_pdf(dist->data,_pars,err);
  forwardError(*err,__LINE__,0);
  return res;
}

void distribution_retrieve(const void* pdist, double* pars, error **err) {
  const distribution *dist;
  
  dist = pdist;
  testErrorRet(dist->retrieve==NULL && dist->n_ded!=0,dist_undef,
	       "Undefined retrieve function for distribution",*err,__LINE__,);
  if (dist->n_ded==0) {
    return;
  }
  dist->retrieve(dist->data,pars,err);
  forwardError(*err,__LINE__,);
}

void free_distribution(distribution **pdist) {
  distribution *dist;
  dist = *pdist;
  
  if (dist->name!=NULL) {
    free(dist->name);
  }
  
  if (dist->free!=NULL && dist->data!=NULL) {
    dist->free(&dist->data);
  }
  if (dist->ndef!=0) {
    free(dist->pars);
    free(dist->def);
  }
  if (dist->dlhandle!=NULL) {
    dlclose(dist->dlhandle);
  }
  free(dist);
  *pdist = NULL;
}

