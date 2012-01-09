/*
 *  clik.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


#include "clik.h"
#include "clik_helper.h"

// ARE YOU STILL READING ?

// YOU HAVE BEEN WARNED !

clik_object* clik_init(char* hdffilepath, error **_err) {
  hid_t file_id,group_id,prior_id;
  herr_t hstat;
  int n_lkl,i_lkl;
  int lmax[6];
  char cur_lkl[100];
  cmblkl **clkl;
  int cli,n_cl;
  zero_bs* zbs;
  distribution *target;
  parname lkl_type;
  
  _dealwitherr;
  
  file_id = H5Fopen( hdffilepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  _testErrorRetVA(file_id<0,hdf5_base,"cannot open  file %s (got %d)",*err,__LINE__,NULL,hdffilepath,file_id);
  
  hstat = H5LTget_attribute_int( file_id, "/clik", "n_lkl_object",  &n_lkl);
  _testErrorRetVA(hstat<0,hdf5_base,"cannot read /clik/n_lkl_object in file %s (got %d)",*err,__LINE__,NULL,hdffilepath,hstat);
  
  hstat = H5LTget_attribute_int( file_id, "/clik", "lmax",  lmax);
  _testErrorRetVA(hstat<0,hdf5_base,"cannot read /clik/lmax in file %s (got %d)",*err,__LINE__,NULL,hdffilepath,hstat);
    
  clkl = malloc_err(sizeof(cmblkl*)*n_lkl,err);
  _forwardError(*err,__LINE__,NULL);
  
  for (i_lkl=0;i_lkl<n_lkl;i_lkl++) {
    sprintf(cur_lkl,"clik/lkl_%d",i_lkl);
    group_id = H5Gopen(file_id, cur_lkl, H5P_DEFAULT );
    _testErrorRetVA(group_id<0,hdf5_base,"cannot read lkl %s in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hdffilepath,hstat);

    clkl[i_lkl] = clik_lklobject_init(group_id,cur_lkl,err);
    _forwardError(*err,__LINE__,NULL);
    
    cmblkl_check_lmax(clkl[i_lkl],lmax,err);
    _forwardError(*err,__LINE__,NULL);
    
    hstat = H5Gclose(group_id);
    _testErrorRetVA(hstat<0,hdf5_base,"cannot close %s in file %s (got %d)",*err,__LINE__,NULL,cur_lkl,hdffilepath,hstat);    
  }
  
  n_cl = 0;
  for(cli=0;cli<6;cli++) {
    n_cl += lmax[cli]+1;
  }
  
  zbs = init_zero_bs(lmax, err);
  _forwardError(*err,__LINE__,NULL);
  
  target = init_multilklbs_distribution(n_cl , clkl,n_lkl,
                                        zbs, &zero_bs_compute, &free_zero_bs, lmax, err);
  _forwardError(*err,__LINE__,NULL);
  
  group_id = H5Gopen(file_id,"clik", H5P_DEFAULT );
  
  hstat = H5Lexists(group_id, "prior", H5P_DEFAULT);
  
  if (hstat==1) {
    int nepar;
    parname *pn;
    int nprior,iprior,j;
    int *lprior;
    char *priorname;
    int nvar;
    double *loc,*var;

    prior_id = H5Gopen(group_id, "prior", H5P_DEFAULT );
    nepar = clik_get_extra_parameter_names(target,&pn,err);
    _forwardError(*err,__LINE__,NULL);
    _testErrorRetVA(nepar==0,hdf5_base,"cannot add a prior without extra parameters",*err,__LINE__,NULL,"");    
    
    nprior = -1;
    priorname = hdf5_char_attarray(prior_id,cur_lkl,"name",&nprior, err);
    _forwardError(*err,__LINE__,NULL);
    nprior = nprior/256;
    _testErrorRetVA(nepar<nprior,hdf5_base,"too many priors ! Expected less than %d got %d",*err,__LINE__,NULL,nepar,nprior);
    
    lprior = malloc_err(sizeof(int)*nprior,err);
    _forwardError(*err,__LINE__,NULL);
    
    for(iprior=0;iprior<nprior;iprior++) {
      lprior[iprior] = -1;
      for (j=0;j<nepar;j++) {
        if (strcmp(pn[j],&(priorname[256*iprior]))==0) {
          lprior[iprior] = n_cl+j;
          break;   
        }
      }
      _testErrorRetVA(lprior[iprior]==-1,hdf5_base,"Unknown extra parameter %s",*err,__LINE__,NULL,priorname[256*iprior]);
    }    
    free(priorname);
    free(pn);
    
    loc = hdf5_double_datarray(prior_id, cur_lkl,"loc",&nprior,err);
    _forwardError(*err,__LINE__,NULL);  
      
    nvar=-1;
    var = hdf5_double_datarray(prior_id, cur_lkl,"var",&nvar,err);
    _forwardError(*err,__LINE__,NULL);  
    
    if (nvar==nprior) {
      target = add_gaussian_prior(target, nprior, lprior, loc, var, err);
      _forwardError(*err,__LINE__,NULL);  
    } else if (nvar==nprior*nprior) {
      target = add_gaussian_prior_2(target, nprior, lprior, loc, var, err);
      _forwardError(*err,__LINE__,NULL);        
    } else {
      _testErrorRetVA(1==1,hdf5_base,"I don't feel well",*err,__LINE__,NULL,"");
    }
    free(loc);
    free(var);
    free(lprior);

    hstat = H5Gclose(prior_id);
    _testErrorRetVA(hstat<0,hdf5_base,"cannot close %s in file %s (got %d)",*err,__LINE__,NULL,"clik/prior",hdffilepath,hstat);    
  
  }


  hstat = H5LTfind_dataset(group_id, "check_param");

  if (hstat==1) {
    int npar;
    double *chkp;
    double res,res2;
    npar = clik_get_extra_parameter_names(target,NULL,err) + n_cl;
    _forwardError(*err,__LINE__,NULL);
    
    chkp = hdf5_double_datarray(group_id, "clik","check_param",&npar,err);
    _forwardError(*err,__LINE__,NULL);
    
    hstat = H5LTread_dataset_double(group_id, "check_value",&res);
    _testErrorRetVA(hstat<0,hdf5_base,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,"check_value","/clik",hstat);
    
    res2 = clik_compute(target,chkp,err);
    _forwardError(*err,__LINE__,NULL);
    
    printf("Checking likelihood '%s' on test data. got %g expected %g (diff %g)\n",hdffilepath,res2,res,res-res2);
    free(chkp);
  }
  hstat = H5Gclose(group_id);
  _testErrorRetVA(hstat<0,hdf5_base,"cannot close %s in file %s (got %d)",*err,__LINE__,NULL,"clik",hdffilepath,hstat);    
  
  hstat = H5Fclose(file_id);
  _testErrorRetVA(hstat<0,hdf5_base,"cannot close ile %s (got %d)",*err,__LINE__,NULL,hdffilepath,hstat);    
  
  return target;
}

void clik_get_has_cl(clik_object *clikid, int has_cl[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  int cli;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,);
  for(cli=0;cli<6;cli++) {
    //fprintf(stderr," %d %d ",cli,lbs->offset_lmax[cli]);
    if (lbs->offset_lmax[cli]!=-1) {
      has_cl[cli]=1;
    } else {
      has_cl[cli]=0;
    }
  }
}

void clik_get_lmax(clik_object *clikid, int lmax[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  zero_bs* zbs;
  int cli;
  _dealwitherr;
  
  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,);
  zbs = lbs->rbs->bs;
  
  for(cli=0;cli<6;cli++) {
    lmax[cli] = zbs->lmax[cli];
  }
}

int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error **_err) {
  parname *pn;
  distribution *target;
  lklbs *lbs;
  int i;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,-1);
  if (names!=NULL) {
    if (lbs->xdim==0) {
      //for now, no extr parameters
      pn = malloc_err(1*sizeof(parname),err);
      _forwardError(*err,__LINE__,-1);
    } else {
      pn = malloc_err(lbs->xdim*sizeof(parname),err);
      _forwardError(*err,__LINE__,-1);
    }
    for(i=0;i<lbs->xdim;i++) {
      sprintf(pn[i],"%s",lbs->xnames[i]);
    }
    *names = pn;  
  }
  return lbs->xdim;
}

int clik_get_extra_parameter_names_by_lkl(clik_object* clikid, int ilkl,parname **names, error **_err) {
  parname *pn;
  distribution *target;
  lklbs *lbs;
  int i;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,-1);
  _testErrorRetVA(ilkl>lbs->nlkl,-11010,"Asked for lkl %d, while there are only %d objects",*err,__LINE__,-1,ilkl,lbs->nlkl);
  
  if (lbs->lkls[ilkl]->xdim==0) {
    //for now, no extr parameters
    pn = malloc_err(1*sizeof(parname),err);
    _forwardError(*err,__LINE__,-1);
  } else {
    pn = malloc_err(lbs->lkls[ilkl]->xdim*sizeof(parname),err);
    _forwardError(*err,__LINE__,-1);
  }
  for(i=0;i<lbs->lkls[ilkl]->xdim;i++) {
    sprintf(pn[i],"%s",lbs->lkls[ilkl]->xnames[i]);
  }
  *names = pn;
  return lbs->lkls[ilkl]->xdim;
}

void clik_cleanup(clik_object** pclikid) {
  free_distribution(pclikid);
}

double clik_compute(clik_object* clikid, double* cl_and_pars,error **_err) {
  double res;
  _dealwitherr;
  
  res = distribution_lkl(clikid, cl_and_pars,err);
  _forwardError(*err,__LINE__,-1);
  return res;
}

void* _clik_dig(clik_object* clikid, error **err) {
  distribution *target;
  target = clikid;
  if (target->log_pdf == &combine_lkl) { 
    // return the first clik likelihood
    int i;
    comb_dist_data* cbd;
    cbd = target->data;
    for (i=0;i<cbd->ndist;i++) {
      if (cbd->dist[i]->log_pdf == &lklbs_lkl) {
        return cbd->dist[i]->data;
      }
    }
  }
  if (target->log_pdf==&lklbs_lkl) {
    return target->data;
  }
  testErrorRet(1==1,-111,"No clik likelihood found",*err,__LINE__,NULL);
}