/*
 *  smica.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "smica.h"

void printMat(double* A, int n, int m) {
  int im,in;
  for(in=0;in<n;in++) {
    for(im=0;im<m-1;im++) {
      fprintf(stderr,"%g , ",A[in*m+im]);
    }
    fprintf(stderr,"%g\n",A[in*m+m-1]);
  }
}
// General funcs

Smica* Smica_init(int nq, double *wq, int m, double *rq_hat, double* rq_0, int nc, SmicaComp **SC,error **err) {
  Smica* smic;
  int isc,iq,info;
  int trois;
  char uplo,diag;
  
  //_DEBUGHERE_("","");
  
  smic = malloc_err(sizeof(Smica),err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  smic->nq = nq;
  smic->m = m;
  trois = 3;
  if(rq_0 == NULL) {
    trois = 2;
  }
  //_DEBUGHERE_("","");
  
  smic->rq_hat = malloc_err(sizeof(double)*(trois*nq+1)*m*m, err);
  forwardError(*err,__LINE__,NULL);
  memcpy(smic->rq_hat,rq_hat,sizeof(double)*m*m*nq);
  
  smic->z_buf = smic->rq_hat + m*m*nq;
  smic->rq = smic->z_buf + m*m;
  //_DEBUGHERE_("","");
  if (rq_0!=NULL) {
    smic->rq_0 = smic->rq + m*m*nq;
    memcpy(smic->rq_0,rq_0,sizeof(double)*m*m*nq);
    //_DEBUGHERE_("","");
  } else {
    smic->rq_0 = NULL;
    //_DEBUGHERE_("","");
  }
  //_DEBUGHERE_("","");
  
  smic->nc = nc;
  smic->SC = malloc_err(sizeof(SmicaComp*)*nc, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  smic->offset_nc = malloc_err(sizeof(int)*nc, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  smic->offset_nc[0] = 0;
  smic->SC[0] = SC[0];
  //_DEBUGHERE_("","");
  for(isc=1;isc<nc;isc++) {
    //_DEBUGHERE_("","");
    smic->offset_nc[isc] = smic->offset_nc[isc-1] + smic->SC[isc-1]->ndim;
    testErrorRetVA(m!=SC[isc]->m,smica_uncomp,"uncompatible number of band in component %d (got %d expected %d)",*err,__LINE__,NULL,isc,SC[isc]->m,m);
    testErrorRetVA(nq!=SC[isc]->nq,smica_uncomp,"uncompatible number of bins in component %d (got %d expected %d)",*err,__LINE__,NULL,isc,SC[isc]->nq,nq);
    
    //_DEBUGHERE_("","");
    smic->SC[isc] = SC[isc];
  }
  //_DEBUGHERE_("","");
  
  smic->wq = malloc_err(sizeof(double)*nq,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  if (wq==NULL) {
    for(iq=0;iq<nq;iq++) {
      //_DEBUGHERE_("","");
      smic->wq[iq] = 1;
    }
  } else {
    //_DEBUGHERE_("","");
    memcpy(smic->wq,wq,sizeof(double)*nq);
  }
  //_DEBUGHERE_("","");
  
  smic->crit = &smica_crit_classic;
  smic->crit_classic_init = 0;
  smic->crit_cor = NULL;
  smic->gvec = NULL;

  smic->eig_buf = NULL;
  smic->eig_lwork = 0;
  smic->eig_nrm = 0;

  smic->lkl_data = NULL;
  smic->lkl_data_free = NULL;
  return smic;
  
}

void free_Smica(void **psmic) {
  Smica *smic;
  int isc;
  
  smic = *psmic;
  
  
  free(smic->wq);
  free(smic->rq_hat);
  free(smic->offset_nc);
  for(isc=0;isc<smic->nc;isc++) {
    if (smic->SC[isc]->names!=NULL) {
      free(smic->SC[isc]->names);
    }
    if (smic->SC[isc]->free!=NULL) {
      smic->SC[isc]->free((void**)&(smic->SC[isc]));          
    }
  }
  free(smic->SC);
    
  if (smic->gvec!=NULL) {
    free(smic->gvec);
  }
  if (smic->quad_mask!=NULL) {
    free(smic->quad_mask);
  }
  if (smic->crit_cor!=NULL) {
    free(smic->crit_cor);
  }
  
  if (smic->eig_buf!=NULL) {
    free(smic->eig_buf);
    free(smic->eig_nrm);
  }
  
  if (smic->lkl_data_free!=NULL) {
    smic->lkl_data_free(&(smic->lkl_data));
  }
  free(smic);
  *psmic=NULL;
}

double Smica_lkl(void* vsmic, double* pars, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  
  smic = vsmic;

  // init rq matrix
  if (smic->rq_0 == NULL) {
    memset(smic->rq,0,sizeof(double)*smic->m*smic->m*smic->nq);
  } else {
    memcpy(smic->rq,smic->rq_0,sizeof(double)*smic->m*smic->m*smic->nq);
  }
  
  // update rq matrix according to each component
  for(isc=0;isc<smic->nc;isc++) {
    char nn[40];
    //_DEBUGHERE_("comp %d update (off %d)",isc,smic->offset_nc[isc]);
    //printMat(smic->rq, smic->m, smic->m);
    //_DEBUGHERE_("%g",*(pars+smic->offset_nc[isc]));
    smic->SC[isc]->update(smic->SC[isc],pars+smic->offset_nc[isc], smic->rq, err);
    forwardError(*err,__LINE__,0);
    //sprintf(nn,"pq_%d.la",isc);
    //write_bin_vector(smic->rq, nn, sizeof(double)*(smic->nq*smic->m*smic->m), err);  
    //forwardError(*err,__LINE__,-1);
  
    //_DEBUGHERE_("comp %d update done",isc);
    //printMat(smic->rq, smic->m, smic->m);
  }
  
  //symetrize matrix
  
  for (iq=0;iq<smic->nq;iq++) {
    double *rq;
    rq = smic->rq+iq*smic->m*smic->m;
    for(i=0;i<smic->m;i++) {
      for (j=0;j<i;j++) {
        rq[j+i*smic->m] = rq[j*smic->m+i];
      }
    }
  }
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);  
  //write_bin_vector(smic->rq_hat, "rqhat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  // ici calculer la vraissemblance a partir de rq et rq_hat
  res = smic->crit(smic,err);
  forwardError(*err,__LINE__,0);

  return res;
}

// modified gaussian criterion
typedef struct {
  int vec_size;
  int *vec_select;
  double *vec_buf,*vec;
  double *sigma_inverse;
} gauss_lkl_data;

double smica_crit_mgauss(void *vsmic, error **err) {
  int iq,im1,im2,iv,m2,m;
  double res;
  Smica *smic;
  char uplo;
  double done,dzero;
  int one,i,j;
  gauss_lkl_data *pld;

  smic = vsmic;
  pld = smic->lkl_data;
  //_DEBUGHERE_("%d %d",smic->nq,smic->m);
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  //write_bin_vector(smic->rq_hat, "rq_hat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  // reorganize data
  //write_bin_vector(smic->quad_mask, "quad_mask.dat", sizeof(int)*(smic->quad_sn), err);   
  m = smic->m;
  m2 = m*m;
  iv = 0;
  for (iv=0;iv<pld->vec_size;iv++) {
    //int civ;
    pld->vec[iv] = smic->rq[pld->vec_select[iv]] - smic->rq_hat[pld->vec_select[iv]];
    //civ = pld->vec_select[iv];
    //_DEBUGHERE_("%d (%d %d %d) %d %g %g %g",iv,civ/m2,(civ-(civ/m2)*m2)/m,(civ-(civ/m2)*m2)%m, pld->vec_select[iv],smic->rq[pld->vec_select[iv]], smic->rq_hat[pld->vec_select[iv]],pld->vec[iv]);
  }
  //_DEBUGHERE_("%d",iv);
  
  one = 1;
  done = 1;
  dzero = 0;
  uplo = 'L';
  //printMat(smic->crit_cor,pld->vec_size,pld->vec_size);
  //_DEBUGHER_("%d",pld->vec_size);
  //write_bin_vector(pld->vec, "gvec.dat", sizeof(double)*(pld->vec_size), err);   
  //write_bin_vector(smic->crit_cor, "crit_cor.dat", sizeof(double)*(pld->vec_size)*(pld->vec_size), err);   
  //_DEBUGHERE_("","");
  //for(i=0;i<pld->vec_size;i++) {
  //  pld->vec_buf[i] = smic->crit_cor[i*pld->vec_size]*pld->vec[0];
  //  for(j=1;j<pld->vec_size;j++) {
  //    pld->vec_buf[i] += smic->crit_cor[i*pld->vec_size+j]*pld->vec[j];
  //  }
  //}
  dsymv(&uplo, &pld->vec_size, &done, pld->sigma_inverse, &pld->vec_size, pld->vec, &one, &dzero, pld->vec_buf, &one);

  //write_bin_vector(pld->vec, "gvecCC.dat", sizeof(double)*(pld->vec_size), err);   
  _DEBUGHERE_("","");
  res = 0;
  for(iq=0;iq<pld->vec_size;iq++) {
    //_DEBUGHERE_("%g %g %g",res,pld->vec[iq],pld->vec_buf[iq])
    res += pld->vec[iq]*pld->vec_buf[iq];
  }
  return -.5*res;  
}

void free_smica_crit_mgauss(void **ppld) {
  gauss_lkl_data *pld;
  pld = *ppld;
  free(pld->vec);
  free(pld->sigma_inverse);
  free(pld);
  _DEBUGHERE_("","");
  *ppld = NULL;
}

void smica_set_crit_mgauss(Smica *smic, int select_size, double *sigma, int *ordering,int* qmin, int* qmax, error **err) {
  int iv,jv,iq,ic;
  int nv,info;
  gauss_lkl_data *pld;
  char uplo;

  pld = malloc_err(sizeof(gauss_lkl_data),err);
  forwardError(*err,__LINE__,);
  _DEBUGHERE_("","");
  
  pld->vec_size = select_size;

  pld->vec_select = malloc_err(sizeof(int)*select_size,err);
  forwardError(*err,__LINE__,);
  
  nv = 0;
  for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
    iv = ordering[ic*2];
    jv = ordering[ic*2+1];
    for (iq=qmin[iv*smic->m+jv];iq<qmax[iv*smic->m+jv];iq++) {
      pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
      nv++;
      testErrorRet(nv>pld->vec_size,-2422,"number of line in covariance matrix smaller than number of selected Cl !",*err,__LINE__,);
    }
  }

  pld->vec = malloc_err(sizeof(double)*pld->vec_size*2,err);
  forwardError(*err,__LINE__,);
  pld->vec_buf = pld->vec + pld->vec_size;
  
  pld->sigma_inverse = malloc_err(sizeof(double)*pld->vec_size*pld->vec_size,err);
  forwardError(*err,__LINE__,);
  memcpy(pld->sigma_inverse,sigma,sizeof(double)*pld->vec_size*pld->vec_size);
  uplo = 'L';
  dpotri(&uplo,&select_size,pld->sigma_inverse,&select_size,&info);
  testErrorRetVA(info!=0,-432432,"cannot inverse covariance with dpotri (info = %d)",*err,__LINE__,,info);

  smic->crit = &smica_crit_mgauss;
  smic->lkl_data = pld;
  smic->lkl_data_free = free_smica_crit_mgauss;
}

//void smica_set_crit_mgauss(Smica *smic, double *crit_cor, int *mask,int *ordering,error **err) {
//  int iv,jv,iq;
//  int nv;
//  gauss_lkl_data *pld;
//
//  pld = malloc_err(sizeof(gauss_lkl_data),err);
//  forwardError(*err,__LINE__,);
//
//  _DEBUGHERE_("","");
//  nv = (smic->nq * smic->m * (smic->m+1))/2;
//  pld->vec_select = malloc_err(sizeof(int)*smic->nq*smic->m*smic->m,err);
//  forwardError(*err,__LINE__,);
//
//  // mask is both in spectra an ell !!!
//  nv = 0;
//  if (ordering == NULL) {
//    for (iv=0;iv<smic->m;iv++) {
//      for (jv=iv;jv<smic->m;jv++) {
//        for (iq=0;iq<smic->nq;iq++) {
//          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
//          if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
//            pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
//            nv++;
//          }
//        }
//      }
//    }
//  } else {
//    int ic;
//    for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
//      iv = ordering[ic*2];
//      jv = ordering[ic*2+1];
//      for (iq=0;iq<smic->nq;iq++) {
//        if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
//          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
//          pld->vec_select[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
//          nv++;
//        }
//      }
//    }
//  }
//
//  pld->vec_size = nv;
//
//  pld->vec = malloc_err(sizeof(double)*pld->vec_size*2,err);
//  forwardError(*err,__LINE__,);
//  pld->vec_buf = pld->vec + pld->vec_size;
//
//  pld->sigma_inverse = malloc_err(sizeof(double)*pld->vec_size*pld->vec_size,err);
//  forwardError(*err,__LINE__,);
//  
//  memcpy(pld->sigma_inverse,crit_cor,sizeof(double)*pld->vec_size*pld->vec_size);
//
//  smic->crit = &smica_crit_mgauss;
//  smic->lkl_data = pld;
//  smic->lkl_data_free = free_smica_crit_mgauss;
//}

// classic criterion

double kld(int n, double* rq_hat, double* rq, error **err) {
  char uplo,trans,diag,side;
  int nn, info,i;
  double *z;
  double done;
  double res;
  
  // suppose que rq_hat et la chol rq_hat
  // rq_hat et rq sont detruits en sortie
  
  //_DEBUGHERE_("","");
  //printMat(rq_hat,n,n);
  //_DEBUGHERE_("","");
  //printMat(rq,n,n);
         
  // chol rq
  uplo = 'L';
  nn = n;
  dpotrf(&uplo,&nn,rq,&nn,&info);
  testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq using dpotrf (%d)",*err,__LINE__,0,info);
  //printMat(rq,nn,nn);
  
  // solve rq z = rq_hat
  side = 'L';
  trans = 'N';
  diag = 'N';
  done = 1;
  dtrsm(&side, &uplo, &trans, &diag, &nn, &nn, &done, rq, &nn, rq_hat, &nn);
  
  z = rq_hat;
  //printMat(z,n,n);
  
  // calcule sum(z^2)
  res = 0;
  //_DEBUGHERE_(": %g",res);
  for(i=0;i<n*n;i++) {
    res += z[i]*z[i];
    //_DEBUGHERE_(": %g",res);
  }
  //_DEBUGHERE_(": %g",res);
  for(i=0;i<n;i++) {
    res -= 2 * log(z[i*n+i]);
    //_DEBUGHERE_(": %g",res);
  }
  //_DEBUGHERE_(": %g",res);
  res -=n;
  //_DEBUGHERE_(": %g",res);
  return 0.5*res;  
}


void smica_set_crit_eig(Smica *smic, double *nrm, error **err) {
  char uplo,jobz;
  int info, m;
  double rq, w,wrk;


  uplo = 'L';
  jobz = 'V';
  info = 0;
  smic->eig_lwork = -1;
  m = smic->m;
  
  dsyev(&jobz, &uplo, &m, &rq, &m, &w, &wrk, &smic->eig_lwork, &info);
  smic->eig_lwork = wrk;

  smic->eig_buf = malloc_err(sizeof(double)*(m+smic->eig_lwork),err);
  forwardError(*err,__LINE__,);

  smic->eig_nrm = malloc_err(sizeof(double)*smic->nq,err);
  forwardError(*err,__LINE__,);

  if (nrm!=NULL) {
    memcpy(smic->eig_nrm,nrm,sizeof(double)*smic->nq);
  } else {
    memset(smic->eig_nrm,0,sizeof(double)*smic->nq);
  }
  
  smic->crit = smica_crit_eig;
  //_DEBUGHERE_("","");
  
}


double smica_crit_eig(void* vsmic, error **err) {
  char jobz,uplo;
  double *wrk, *w;
  int info;
  double res,pes;
  int i,j,k,m,lwork,iq;
  Smica *smic;
  double jes,les,dg;
  double *rq,*rq_hat;
  smic = vsmic;

  m = smic->m;
  // assumes rq_hat is rq_hat !
  // nrm is -.5 (log(det(rq_hat)) + ndim) or any arbitrary number
  // buf is big enough to hold wrk and w
  // compute the eigenvalues and modes of rq

  uplo = 'L';
  jobz = 'V';
  info = 0;
  w = smic->eig_buf;
  wrk = w + m;
  lwork = smic->eig_lwork;
  res = 0;
  
  for (iq=0;iq<smic->nq;iq++) {
    rq = &(smic->rq[iq*m*m]);
    rq_hat = &(smic->rq_hat[iq*m*m]);

    dsyev(&jobz, &uplo, &m, rq, &m, w, wrk, &lwork, &info);
    testErrorRetVA(info !=0,-1234321,"dsyev failed (%d)",*err,__LINE__,0,info);
  
    les = 0;

    // ici peut etre un peu de regul du resultat, un jour
    for(i=0;i<m;i++) {
      les += log(w[i]); 
      pes = 0;
      dg=0;
      for(j=0;j<m;j++) {
        jes = 0;
        for(k=0;k<m;k++) {
          jes += rq_hat[j*m+k] * rq[k+i*m];
        }
        dg += rq[i+j*m] * rq[i+j*m];
        pes += jes * rq[j+i*m];
      }
      les += pes/w[i];
    }

    //_DEBUGHERE_("%g",les*.5 - smic->eig_nrm[iq])
    
    res += (les*.5 - smic->eig_nrm[iq]) * smic->wq[iq];
}
  
  return -res;
}


double smica_crit_classic(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;


  smic = vsmic;
  m = smic->m;
  nq = smic->nq;
  if (smic->crit_classic_init==0) {
    for(iq = 0; iq < nq;iq++) { //precompute chol decomposed rq_hat
      int mx,my,info;
      double *rql;
      char uplo;
      //_DEBUGHERE_("","");
      
      rql = smic->rq_hat + iq*m*m;
      // chol
      uplo = 'L';
      dpotrf(&uplo,&m,rql,&m,&info);
      testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq_hat using dpotrf (%d)",*err,__LINE__,0,info);
      //_DEBUGHERE_("","");
      
      // fill the U part with 0 (beware, f90!)
      for(mx=0;mx<m;mx++) {
        for(my=mx+1;my<m;my++) {
          rql[my*m+mx] = 0;
        }
      }
      //_DEBUGHERE_("","");  
    }
    smic->crit_classic_init=1;
  }

  res = 0;
  for(iq=0;iq<nq;iq++) {
    double kdd;
    //_DEBUGHERE_("iq %d -> %g",iq,res);
    memcpy(smic->z_buf,smic->rq_hat+m*m*iq,m*m*sizeof(double));
    kdd = kld(m,smic->z_buf,smic->rq+m*m*iq,err);
    //_DEBUGHERE_("%g",kdd);
    res += smic->wq[iq] * kdd;
    //_DEBUGHERE_(" -> %g (%g)",res,smic->wq[iq]);
    forwardError(*err,__LINE__,0);
  }
  //_DEBUGHERE_(" --> %g (%g)",res);  
  return -res; 
}

// quad criterion

void smica_set_crit_quad(Smica *smic, double *fid,int *mask,error **err) {
  int iq,m,nq;
  char uplo,trans,diag,side;
  int info;
  double *fdq;
  int sz;
  int i,j,sn;

  //_DEBUGHERE_("","");
  smic->crit = smica_crit_quad;
  smic->quad_sn = smic->m;
  smic->quad_mask = NULL;
  if (mask!=NULL) {
    smic->quad_mask = malloc_err(sizeof(int)*smic->m*smic->m,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->quad_mask,mask,sizeof(int)*smic->m*smic->m);
    sn = 0;
    
    for(i=0;i<smic->m;i++) {
      for(j=i;j<smic->m;j++) {
        sn += mask[i*smic->m+j];
      }
    }
    smic->quad_sn = sn;
    smic->eig_buf = malloc_err(sizeof(double)*sn*sn*smic->nq,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->eig_buf,fid,sizeof(double)*sn*sn*smic->nq);
    smic->crit = smica_crit_quad_mask;
    
    return;
  }
  if (fid!=NULL) {
    //_DEBUGHERE_("","");
    smic->crit = smica_crit_quadfid;
    smic->eig_buf = malloc_err(sizeof(double)*smic->m*smic->m*smic->nq,err);
    forwardError(*err,_LINE__,);
    memcpy(smic->eig_buf,fid,sizeof(double)*smic->m*smic->m*smic->nq);
    uplo = 'L';
    m = smic->m;
    sz = m*m;

    for(iq = 0;iq<smic->nq;iq++) {
      fdq = smic->eig_buf + iq*sz;
      dpotrf(&uplo,&m,fdq,&m,&info);
      testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq_hat using dpotrf (%d)",*err,__LINE__,,info);   
    }
  }
}

double smica_crit_quad_mask(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq,*c_rq_fid;
  char uplo,trans,diag,side;
  int info;
  int i,j,p;
  double *vecq;
  double done,dzero;
  int sn;
  double kdd;

  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
  done = 1;
  dzero = 1;
  sn = smic->quad_sn;

  res = 0;
      
  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->rq_hat + sz*iq;
    c_rq = smic->rq + sz*iq;
    c_rq_fid = smic->eig_buf + sn*sn*iq;
    vecq = smic->z_buf;

    //printMat(c_rq,m,m);
    // c_rq = c_rq - c_rq_hat  
    daxpy(&sz, &dminusone,c_rq_hat,&one, c_rq, &one);
    
    p = 0;
    for(i=0;i<m;i++) {
      for(j=i;j<m;j++) {
        
        if (smic->quad_mask[i*m+j]==0) {
          continue;
        }
        if (i==j) {
          vecq[p] = c_rq[i*m+j];
        } else {
          vecq[p] = c_rq[i*m+j];
        }
        
        p++;
      }
    }
    
    
    dsymv(&uplo, &sn, &done, c_rq_fid, &sn, vecq, &one, &dzero, c_rq, &one);
    
        
    kdd = 0;
    for(i=0;i<sn;i++) {
      kdd += c_rq[i]*vecq[i];
    }
    
    kdd=kdd/4.;
    

    res += smic->wq[iq] * kdd;
  }

  return -res;  
}

double smica_crit_quadfid(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq,*c_rq_fid;
  char uplo,trans,diag,side;
  int info;
  int i,j;

  
  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
    
  res = 0;

  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->rq_hat + sz*iq;
    c_rq = smic->rq + sz*iq;
    c_rq_fid = smic->eig_buf + sz*iq;
    
    //printMat(c_rq,m,m);
    // c_rq = c_rq - c_rq_hat  
    daxpy(&sz, &dminusone,c_rq_hat,&one, c_rq, &one);
    

    dpotrs(&uplo, &m, &m, c_rq_fid, &m, c_rq, &m, &info);
    testErrorRetVA(info!=0,lowly_chol,"Could not solve rq using dpotrs (%d)",*err,__LINE__,0,info);

    kdd = 0;
    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        kdd += c_rq[i*m+j]*c_rq[j*m+i];
      }
    }
    kdd=kdd/4.;
    
    res += smic->wq[iq] * kdd;
  }

  return -res;  
}


double smica_crit_quad(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;
  int one,sz;
  double dminusone;
  double *c_rq_hat,*c_rq;
  char uplo,trans,diag,side;
  int info;
  int i,j;

  smic = vsmic;
  m = smic->m;
  nq = smic->nq;

  one = 1;
  sz = m*m;
  dminusone = -1;
  uplo = 'L';
    
  res = 0;

  for(iq=0;iq<nq;iq++) {
    double kdd;
    c_rq_hat = smic->z_buf;
    c_rq = smic->rq + sz*iq;
    memcpy(c_rq_hat,smic->rq_hat+sz*iq,m*m*sizeof(double));
    
    // c_rq_hat = c_rq_hat - c_rq  
    daxpy(&sz, &dminusone,c_rq,&one, c_rq_hat, &one);

    dpotrf(&uplo,&m,c_rq,&m,&info);
    testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq using dpotrf (%d)",*err,__LINE__,0,info);
  
    dpotrs(&uplo, &m, &m, c_rq, &m, c_rq_hat, &m, &info);
    testErrorRetVA(info!=0,lowly_chol,"Could not solve rq using dpotrs (%d)",*err,__LINE__,0,info);

    kdd = 0;
    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        kdd += c_rq_hat[i*m+j]*c_rq_hat[j*m+i];
      }
    }
    kdd=kdd/4.;

    res += smic->wq[iq] * kdd;
  }

  return -res;  
}

//gaussian approx criterion

double smica_crit_gauss(void *vsmic, error **err) {
  int iq,im1,im2,iv,m2,m;
  double res;
  Smica *smic;
  char uplo;
  double done,dzero;
  int one,i,j;

  smic = vsmic;
  //_DEBUGHERE_("%d %d",smic->nq,smic->m);
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  //write_bin_vector(smic->rq_hat, "rq_hat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  // reorganize data
  //write_bin_vector(smic->quad_mask, "quad_mask.dat", sizeof(int)*(smic->quad_sn), err);   
  m = smic->m;
  m2 = m*m;
  iv = 0;
  for (iv=0;iv<smic->quad_sn;iv++) {
    int civ;
    smic->gvec[iv] = smic->rq[smic->quad_mask[iv]] - smic->rq_hat[smic->quad_mask[iv]];
    civ = smic->quad_mask[iv];
    //_DEBUGHERE_("%d (%d %d %d) %d %g %g %g",iv,civ/m2,(civ-(civ/m2)*m2)/m,(civ-(civ/m2)*m2)%m, smic->quad_mask[iv],smic->rq[smic->quad_mask[iv]], smic->rq_hat[smic->quad_mask[iv]],smic->gvec[iv]);
  }
  //_DEBUGHERE_("%d",iv);
  
  one = 1;
  done = 1;
  dzero = 0;
  uplo = 'L';
  //printMat(smic->crit_cor,smic->quad_sn,smic->quad_sn);
  //_DEBUGHER_("%d",smic->quad_sn);
  //write_bin_vector(smic->gvec, "gvec.dat", sizeof(double)*(smic->quad_sn), err);   
  //write_bin_vector(smic->crit_cor, "crit_cor.dat", sizeof(double)*(smic->quad_sn)*(smic->quad_sn), err);   
  //_DEBUGHERE_("","");
  //for(i=0;i<smic->quad_sn;i++) {
  //  smic->gvec[smic->quad_sn+i] = smic->crit_cor[i*smic->quad_sn]*smic->gvec[0];
  //  for(j=1;j<smic->quad_sn;j++) {
  //    smic->gvec[smic->quad_sn+i] += smic->crit_cor[i*smic->quad_sn+j]*smic->gvec[j];
  //  }
  //}
  dsymv(&uplo, &smic->quad_sn, &done, smic->crit_cor, &smic->quad_sn, smic->gvec, &one, &dzero, smic->gvec+smic->quad_sn, &one);

  //write_bin_vector(smic->gvec, "gvecCC.dat", sizeof(double)*(smic->quad_sn), err);   
  
  res = 0;
  for(iq=0;iq<iv;iq++) {
    //_DEBUGHERE_("%g %g %g",res,smic->gvec[iq],smic->gvec[iq+iv])
    res += smic->gvec[iq]*smic->gvec[iq+iv];
  }
  return -.5*res;  
}

void smica_set_crit_gauss(Smica *smic, double *crit_cor, int *mask,int *ordering,error **err) {
  int iv,jv,iq;
  int nv;


  nv = (smic->nq * smic->m * (smic->m+1))/2;
  smic->quad_mask = malloc_err(sizeof(int)*smic->nq*smic->m*smic->m,err);
  forwardError(*err,__LINE__,);

  // mask is both in spectra an ell !!!
  nv = 0;
  if (ordering == NULL) {
    for (iv=0;iv<smic->m;iv++) {
      for (jv=iv;jv<smic->m;jv++) {
        for (iq=0;iq<smic->nq;iq++) {
          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
          if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
            smic->quad_mask[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
            nv++;
          }
        }
      }
    }
  } else {
    int ic;
    for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
      iv = ordering[ic*2];
      jv = ordering[ic*2+1];
      for (iq=0;iq<smic->nq;iq++) {
        if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
          smic->quad_mask[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
          nv++;
        }
      }
    }
  }

  smic->gvec = malloc_err(sizeof(double)*nv*2,err);
  forwardError(*err,__LINE__,);
  smic->quad_sn = nv;

  smic->crit_cor = malloc_err(sizeof(double)*nv*nv,err);
  forwardError(*err,__LINE__,);
  
  memcpy(smic->crit_cor,crit_cor,sizeof(double)*nv*nv);

  smic->crit = &smica_crit_gauss;
}


// Simple components

// constant

SmicaComp* comp_cst_init(int nq, int m, double *rq_0, error **err) {
  SmicaComp *SC;
  double *data;

  data = malloc_err(sizeof(double)*m*m*nq,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(data,rq_0,sizeof(double)*m*m*nq);

  SC = alloc_SC(0,nq,m,data,&comp_cst_update,&free_comp_cst,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}

void comp_cst_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  double * rq_0;
  int sz,one;
  double done;

  SC = data;
  rq_0 = SC->data;
  sz = SC->m*SC->m*SC->nq;
  done = 1;
  one = 1;

  daxpy(&sz,&done,rq_0,&one, rq, &one);
}

void free_comp_cst(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  free(SC->data);
  free(SC);
  *data = NULL;
}

// 1D

void comp_1D_AAt(int m, double *A, double *AAt, error **err) {
  int mm,one;
  char uplo;
  double done;
    
  memset(AAt,0,m*m*sizeof(double));
  
  mm = m;
  uplo = 'L';
  done = 1;
  one = 1;
  // -> AAt = 1 * A * A' + AAt
  //printMat(A, m, 1);
  dsyr(&uplo, &mm, &done, A, &one, AAt, &mm);  
  //printMat(AAt, m, m);
}





SmicaComp* alloc_SC(int ndim,int nq,int m,void* data, update_rq* update, posterior_log_free* pfree, error **err) {
  SmicaComp* SC;
  
  SC = malloc_err(sizeof(SmicaComp), err);
  forwardError(*err,__LINE__,NULL);
  SC->m = m;
  SC->nq = nq;
  SC->ndim = ndim;
  SC->update = update;
  SC->data = data;
  SC->free = pfree;
  SC->names=NULL;
  SC_set_compname(SC,"UNK");
  return SC;
}

void SC_set_compname(SmicaComp *SC, char *name) {
  sprintf(SC->comp_name,"%s",name);
}


void SC_setnames(SmicaComp *SC, char** names, error **err) {
  int i;
  if (SC->names!=NULL) {
    free(SC->names);
  }
  if (SC->ndim!=0) {
   SC->names = malloc_err(sizeof(_smicanames)*SC->ndim,err);
   forwardError(*err,__LINE__,);  
  } else{
    SC->names = malloc_err(sizeof(_smicanames)*1,err);
   forwardError(*err,__LINE__,);  
  }
  for(i=0;i<SC->ndim;i++) {
    sprintf(SC->names[i],"%s",names[i]);
  }
}


SmicaComp* comp_1D_init(int nq, int m, double *A, error **err) {
  SmicaComp *SC;
  char uplo;
  double done;
  int one;
  SC_1D_data *data;
  int ndim;
  
  data = malloc_err(sizeof(SC_1D_data), err);
  forwardError(*err,__LINE__,NULL);
  
  data->AAt = malloc_err(sizeof(double)*m*m, err);
  forwardError(*err,__LINE__,NULL);
  
  if (A!=NULL) {
    data->Acst=1;
    ndim = nq;

    comp_1D_AAt(m, A,data->AAt, err);    
    forwardError(*err,__LINE__,NULL);
    
  } else {
    data->Acst=0;
    ndim = nq + m;    
  }
  
  SC = alloc_SC(ndim,nq,m,data,&comp_1D_update,&free_comp_1D,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}




void free_comp_1D(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  //_DEBUGHERE_("","");
  free(((SC_1D_data*) SC->data)->AAt);
  //_DEBUGHERE_("","");
  free(SC->data);
  //_DEBUGHERE_("","");
  free(SC);
  //_DEBUGHERE_("","");
  *data = NULL;
  //_DEBUGHERE_("","");
}




void comp_1D_update(void* data,double* locpars, double* rq, error **err) {
  int iq,one;
  SmicaComp *SC;
  double *AAt;
  int m,m2;
  double *mpars;
  SC_1D_data *SCdat;
  
  SC = data;
  m = SC->m;
  SCdat = SC->data;
  AAt = SCdat->AAt;
  m2 = m*m;
  one = 1;
  
  mpars = locpars;
  if (SCdat->Acst==0) {
    double *A;
    A = locpars;
    mpars = locpars + m;
    comp_1D_AAt(m,A,AAt, err);    
    forwardError(*err,__LINE__,);
  }
  
  //printMat(AAt, m, m);

  for(iq=0;iq<SC->ndim;iq++) {
    //printMat(rq+m2*iq, m, m);
    daxpy(&m2,mpars + iq,AAt,&one, rq+m2*iq, &one);
    //printMat(rq+m2*iq, m, m);
    // -> rq[iq*m2] = locpars[iq] * AAt + rq[iq*m2]
  }
}



// nD

SmicaComp* comp_nD_init(int nq, int m, int nd, double *A, error **err) {
  SmicaComp *SC;
  char uplo;
  double done;
  int one;
  int ndim;
  void* data;
  
  
  //_DEBUGHERE_("","");
  data = malloc_err(sizeof(SC_nD_data), err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->nd = nd;

  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->A = malloc_err(sizeof(double)*(m*nd*2+nd*nd), err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  if (A!=NULL) {
    ndim = nq*(nd*(nd+1))/2;
    ((SC_nD_data*) data)->Acst=1;
    //_DEBUGHERE_("","");
    
    memcpy(((SC_nD_data*) data)->A,A,m*nd*sizeof(double));    
  } else {
    ndim = m*nd+nq*(nd*(nd+1))/2;
    ((SC_nD_data*) data)->Acst=0;
    //_DEBUGHERE_("","");
    
  }
  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->Ab = ((SC_nD_data*) data)->A + m*nd;
  ((SC_nD_data*) data)->P = ((SC_nD_data*) data)->Ab + m*nd;
  //_DEBUGHERE_("","");
  
  // beware A is C oriented
  //_DEBUGHERE_("","");
  SC = alloc_SC(ndim,nq,m,data,&comp_nD_update,&free_comp_nD,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  return SC;
}




void free_comp_nD(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  //_DEBUGHERE_("","");
  free(((SC_nD_data*) SC->data)->A);
  //_DEBUGHERE_("","");
  free(SC->data);
  //_DEBUGHERE_("","");
  free(SC);
  //_DEBUGHERE_("","");
  *data = NULL;
}




void comp_nD_update(void* data,double* locpars, double* rq, error **err) {
  int iq,one;
  SmicaComp *SC;
  double *A,*Ab,*P,*Ppack;
  int m,m2,nd,nq;
  double done,dzero;
  double *mpars;
  char transa,transb,side,uplo;
  
  // locpar is a (C oriented) q*nd*nd U Triangular matrix
  
  SC = data;
  m = SC->m;
  nq = SC->nq;
  
  nd = ((SC_nD_data*) SC->data)->nd;
  Ab = ((SC_nD_data*) SC->data)->Ab;
  P = ((SC_nD_data*) SC->data)->P;
  if (((SC_nD_data*) SC->data)->Acst==1) {
    A = ((SC_nD_data*) SC->data)->A;
    mpars = locpars;
  } else {
    A = locpars;
    mpars = locpars + m*nd;
  }
  
  m2 = m*m;
  one = 1;
  
  done = 1;
  dzero = 0;
  transb = 'N'; // because A is C ordered
  
  for(iq=0;iq<nq;iq++) {
    int ii;
    int ix,iy;
    
    // unpack locpar[q]
    Ppack = locpars + iq*((nd*(nd+1))/2);
    ii = 0;
    for(ix=0;ix<nd;ix++) {
      for(iy=ix;iy<nd;iy++) {
        P[ix*nd+iy] = Ppack[ii];
        //P[iy*nd+ix] = Ppack[ii];
        //P[ix*nd+iy] = 0;
        //P[iy*nd+ix] = 0;
        ii++;
      }
      //P[ix*nd+ix] = 1;
    }
    //printMat(P,2,2);

    //_DEBUGHERE_("P %d",iq);
    //printMat(P, nd, nd);
    //_DEBUGHERE_("A","");
    //printMat(A, m, nd);
    //_DEBUGHERE_("P.A'","");
    //printMat(Ab, m, nd);    
    transa = 'N';
    side = 'L';
    uplo = 'L';
    // Ab = P.A' (Ab is fortran ordered while A is C ordered)
    dsymm(&side, &uplo, &nd, &m, &done, P, &nd, A, &nd, &dzero, Ab, &nd);
    //dgemm(&transa, &transb, &nd, &m, &nd, &done, P, &nd, A, &nd, &dzero, Ab, &nd);

    //_DEBUGHERE_("P %d",iq);
    //printMat(P, nd, nd);
    //_DEBUGHERE_("A","");
    //printMat(A, m, nd);    
    //_DEBUGHERE_("P.A'","");
    //printMat(Ab, m, nd);    
    
    // Rq += A.Ab (=A.P.A')
    transa = 'T'; // because A is c ordered

    //_DEBUGHERE_("rq a","");
    //printMat(rq+m2*iq, m, m);    
    dgemm(&transa, &transb, &m, &m, &nd, &done, A, &nd, Ab, &nd, &done, rq+m2*iq, &m);
    //_DEBUGHERE_("rq b","");
    //printMat(rq+m2*iq, m, m);    
    //printMat(rq+m2*iq,2,2);
  }
}



//diagonal

SmicaComp* comp_diag_init(int nq, int m, error **err) {
  SmicaComp* SC;
  
  SC = alloc_SC(nq*m,nq,m,NULL,&comp_diag_update,NULL,err);
  forwardError(*err,__LINE__,NULL);
  return SC;
}




void comp_diag_update(void* data,double* locpars, double* rq, error **err) {
  int iq,im,m2,ii,m;
  SmicaComp* SC;
  
  SC = data;
  m = SC->m;
  
  m2=m*m;
  ii=0;
  //_DEBUGHERE_("pars","");
  //printMat(locpars, SC->nq, m);    
  for(iq=0;iq<SC->nq;iq++) {
    //_DEBUGHERE_("rq a","");
    //printMat(rq+m2*iq, m, m);    
    for(im=0;im<m;im++) {
      rq[iq*m2+im*m+im] += locpars[ii];
      ii++;
    }
    //_DEBUGHERE_("rq b","");
    //printMat(rq+m2*iq, m, m);    
  }
}

SmicaComp* amp_diag_init(int nq, int m, double* tmpl, error **err) {
  SmicaComp *SC;
  double *itmpl;

  itmpl = malloc_err(sizeof(double)*nq*m,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(itmpl,tmpl,sizeof(double)*nq*m);

  SC = alloc_SC(m,nq,m,itmpl,&amp_diag_update,&amp_diag_free,err);
  forwardError(*err,__LINE__,NULL);
  
  return SC;
}

void amp_diag_update(void* data,double* locpars, double* rq, error **err) {
  int iq,im,m2,ii,m,q2;
  SmicaComp* SC;
  double *tmpl;

  SC = data;
  m = SC->m;
  tmpl = SC->data;
  m2 = m*m;
  q2 = SC->nq*SC->nq;

  ii=0;
  //_DEBUGHERE_("pars","");
  //printMat(locpars, SC->nq, m);
  for(im=0;im<SC->m;im++) {
    double pp;
    pp = locpars[im];
    for(iq=0;iq<SC->nq;iq++) {
      rq[iq*m2+im*m+im] += pp*tmpl[im*q2+iq];
    }
  }     
}

void amp_diag_free(void** data) {
  SmicaComp *SC;
  SC = *data;
  free(SC->data);
  free(SC);
  *data = NULL;
}

// CMB

SmicaComp * comp_CMB_init(int nbins, int mt,int mp, int *has_cl, double* Acprs, error **err) {
  double *A;
  SC_CMB_data* data;
  int mtot;
  int trois,im,i,six;
  SmicaComp *SC;
  
  data = malloc_err(sizeof(SC_CMB_data), err);
  forwardError(*err,__LINE__,NULL);
  
  trois = has_cl[0]+has_cl[1]+has_cl[2];
  testErrorRet(trois==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  six = trois + has_cl[3]+has_cl[4]+has_cl[5];
  
  mtot = mt*has_cl[0]+(has_cl[1]+has_cl[2])*mp;
  testErrorRet(mtot==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  
  testErrorRet(mt!=0 && has_cl[0]==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(mt==0 && has_cl[0]!=0,smica_uncomp,"mismatch",*err,__LINE__,NULL);

  testErrorRet(mp!=0 && has_cl[1]==0 && has_cl[2]==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(mp==0 && (has_cl[1]!=0 || has_cl[2]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);

  testErrorRet(has_cl[0]==0 && (has_cl[3]!=0 || has_cl[4]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(has_cl[1]==0 && (has_cl[3]!=0 || has_cl[5]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(has_cl[2]==0 && (has_cl[4]!=0 || has_cl[5]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);

  if (trois==1) {
    // cas particulier 1D
    data->locpars=NULL;

    data->SCnD = comp_1D_init(nbins,mtot,Acprs,err);
    forwardError(*err,__LINE__,NULL);
  
    SC = alloc_SC(nbins,nbins,mtot,data,&comp_CMB_update,&free_comp_CMB,err);
    forwardError(*err,__LINE__,NULL);
    return SC;
    
  }
  
  A = malloc_err(sizeof(double)*mtot*trois, err);
  forwardError(*err,__LINE__,NULL);
  memset(A,0,sizeof(double)*mtot*trois);
  // fill T
  if(has_cl[0]==1) {
    for(im=0;im<mt;im++) {
      A[im*trois] = Acprs[im];
    }
  }
  // fill E and B
  for(im=0;im<mp;im++) {
    if (has_cl[1]==1) {
      A[(im+mt*has_cl[0])*trois+has_cl[0]] = Acprs[im+mt];      
    }
    if (has_cl[2]==1) {
      A[(im+mt*has_cl[0]+has_cl[1]*mp)*trois+has_cl[0]+has_cl[1]] = Acprs[im+mt];      
    }
  }
  
  data->locpars = malloc_err(sizeof(double)*nbins*(trois*(trois+1))/2,err);
  forwardError(*err,__LINE__,NULL);
  
  data->SCnD = comp_nD_init(nbins,mtot,trois,A,err);
  forwardError(*err,__LINE__,NULL);

  free(A);

  for(i=0;i<6;i++) {
    data->has_cl[i] = has_cl[i];
    data->jmp_cl[i] = -1;
  }
  
  if (trois==3) {
    data->jmp_cl[0] = 0;
    data->jmp_cl[1] = 3;
    data->jmp_cl[2] = 5;
    if (has_cl[3]==1) {
      data->jmp_cl[3] = 1;
    }
    if (has_cl[4]==1) {
      data->jmp_cl[4] = 2;
    }
    if (has_cl[5]==1) {
      data->jmp_cl[5] = 4;
    }
  }
  
  if (trois==2) {
    if (has_cl[0]==1) {
      data->jmp_cl[0] = 0;
      if (has_cl[1]==1) {
        data->jmp_cl[1] = 2;
        if (has_cl[3]==1) {
          data->jmp_cl[3] = 1;  
        }
      } else {
        data->jmp_cl[2] = 2;
        if (has_cl[4]==1) {
          data->jmp_cl[4] = 1;
        }
      }
    } else {
      data->jmp_cl[1] = 0;
      data->jmp_cl[2] = 2;
      if (has_cl[5]==1) {
        data->jmp_cl[5] = 1;
      }
    }
  }
  data->trois = trois;
  
  SC = alloc_SC(nbins*six,nbins,mtot,data,&comp_CMB_update,&free_comp_CMB,err);
  forwardError(*err,__LINE__,NULL);
  return SC;
  
}




void comp_CMB_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_CMB_data *SCd;
  int td,ic,im,i,iq;
  
  SC = data;
  SCd = SC->data;
  
  if (SCd->locpars==NULL) {
    // 1d case
    comp_1D_update(SCd->SCnD,locpars,rq,err);
    forwardError(*err,__LINE__,);
    return;
  } 
  
  td  = (SCd->trois*(SCd->trois+1))/2;
  
  // nd il faut que je reoriente mon vecteur
  i=0;
  for(ic=0;ic<6;ic++) {

    if(SCd->jmp_cl[ic]!=-1) {
      for(iq=0;iq<SC->nq;iq++) {
        //_DEBUGHERE_("%d %d %d %g %d",ic,iq,i, locpars[i],iq*td+SCd->jmp_cl[ic]);
        SCd->locpars[iq*td+SCd->jmp_cl[ic]] = locpars[i];
        i+=SCd->has_cl[ic];
      }
    }
  }
  comp_nD_update(SCd->SCnD,SCd->locpars,rq,err);
  forwardError(*err,__LINE__,);
  return;
  
}




void free_comp_CMB(void** data) {
  SmicaComp *SC;
  SC_CMB_data *SCd;
  
  SC = *data;
  SCd = SC->data;
  SCd->SCnD->free((void**)&SCd->SCnD);
  if(SCd->locpars!=NULL) {
    free(SCd->locpars);
  }
  free(SCd);
  free(SC);
  *data = NULL;
}


SC_gcal* comp_gcal_gen_init(int q, int m, int* ngcal, double *gcaltpl, int nell, double *bins, error **err) {
  SC_gcal *gc;
  int ntpl,im;
  
  gc = malloc_err(sizeof(SC_gcal),err);
  forwardError(*err,__LINE__,NULL);
  
  gc->ngcal = malloc_err(sizeof(int)*m,err);
  forwardError(*err,__LINE__,NULL);
  
  ntpl=0;
  for(im=0;im<m;im++) {
    ntpl+=ngcal[im];
    gc->ngcal[im]=ngcal[im];
  }
  gc->ntpl = ntpl;
  
  gc->gcaltpl = malloc_err(sizeof(double)*q*ntpl,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*q*ntpl);

  if (nell>0) {
    int iq,il;
    gc->gcaltpl = malloc_err(sizeof(double)*nell*ntpl,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*nell*ntpl);
    gc->nell = nell;
    gc->tpl = malloc_err(sizeof(double)*nell*m,err);
    forwardError(*err,__LINE__,NULL); 
    gc->bins = malloc_err(sizeof(double)*(nell*q),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->bins,bins,sizeof(double)*q*nell);
    gc->tpll = malloc_err(sizeof(double)*q,err);
    forwardError(*err,__LINE__,NULL);
    for (iq=0;iq<q;iq++) {
      gc->tpll[iq] = 0;
      for(il=0;il<nell;il++) {
        gc->tpll[iq] += bins[iq*nell+il]*bins[iq*nell+il];
      }
    }    
  } else {
    gc->nell =0;
    gc->gcaltpl = malloc_err(sizeof(double)*q*ntpl,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(gc->gcaltpl,gcaltpl,sizeof(double)*q*ntpl);
    gc->tpl = malloc_err(sizeof(double)*q*m,err);
    forwardError(*err,__LINE__,NULL);
    gc->tpll=NULL;
    gc->bins = NULL;
  }
  return gc;
}

void comp_gcal_apply_update(int m,int nq,double* rq,double *tpl,int nell,double *tpll, double *bins) {
  int iq,im1,im2,il;

  if (nell>0) {
    for (iq=0;iq<nq;iq++) {
      for(im1=0;im1<m;im1++) {  
        for(im2=0;im2<m;im2++) {
          double cr;
          cr=0;
          for(il=0;il<nell;il++) {
            cr += bins[iq*nell+il] * bins[iq*nell+il] * tpl[im1*nq + il] * tpl[im2*nq + il];
          }
          rq[iq*m*m+im1*m+im2] *= cr/tpll[iq];
        }
      }
    }  
  } else {
    for (iq=0;iq<nq;iq++) {
      for(im1=0;im1<m;im1++) {  
        for(im2=0;im2<m;im2++) {
          rq[iq*m*m+im1*m+im2] *= tpl[im1*nq + iq] * tpl[im2*nq + iq];
        }
      }
    }
  }
}


SmicaComp* comp_gcal_lin_init(int q,int m, int *ngcal, double* gcaltpl, int nell, double *bins,error **err) {
  SC_gcal *gc;
  SmicaComp *SC;

  gc = comp_gcal_gen_init(q,m,ngcal,gcaltpl,nell,bins,err);
  forwardError(*err,__LINE__,NULL);

  SC = alloc_SC(gc->ntpl,q,m,gc, &comp_gcal_lin_update, &comp_gcal_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}

void comp_gcal_lin_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal *gc;
  int im1,im2, iq,ig;
  int m2,m,nq,ipars;

  SC = data;
  gc = SC->data;
  m2 = SC->m*SC->m;
  m = SC->m;
  nq = SC->nq;
  ipars = 0;

  if (gc->nell>0) {
    nq = gc->nell;
  }
  ipars = 0;
  // first build the template 
  for(im1=0;im1<m;im1++) {
    if (gc->ngcal[im1]==0) {
      for(iq=0;iq<nq;iq++) {
          gc->tpl[im1*nq + iq] = 1;
      }
      continue;
    }
    for(iq=0;iq<nq;iq++) {
      gc->tpl[im1*nq + iq] = locpars[ipars] * gc->gcaltpl[ipars*nq+iq];
    }
    ipars++;
    for(ig=1;ig<gc->ngcal[im1];ig++) {
      for(iq=0;iq<nq;iq++) {
        gc->tpl[im1*nq + iq] += locpars[ipars] * gc->gcaltpl[ipars*nq+iq];
       
      }
      ipars++;
    }
  }

  comp_gcal_apply_update(m,SC->nq,rq,gc->tpl,gc->nell,gc->tpll,gc->bins);
}

SmicaComp* comp_gcal_log_init(int q,int m, int* ngcal,  double* gcaltpl,int nell, double *bins,error **err) {
  SC_gcal *gc;
  SmicaComp *SC;

  gc = comp_gcal_gen_init(q,m,ngcal,gcaltpl,nell,bins,err);
  forwardError(*err,__LINE__,NULL);

  SC = alloc_SC(gc->ntpl,q,m,gc, &comp_gcal_log_update, &comp_gcal_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}

void comp_gcal_log_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal *gc;
  int im1,im2, iq,ig;
  int m2,m,nq,ipars;

  SC = data;
  gc = SC->data;
  m2 = SC->m*SC->m;
  m = SC->m;
  nq = SC->nq;
  ipars = 0;

  if (gc->nell>0) {
    nq = gc->nell;
  }
  ipars = 0;
  // first build the template 
  for(im1=0;im1<m;im1++) {
    if (gc->ngcal[im1]==0) {
      for(iq=0;iq<nq;iq++) {
          gc->tpl[im1*nq + iq] = 1;
      }
      continue;
    }
    for(iq=0;iq<nq;iq++) {
      gc->tpl[im1*nq + iq] = exp(locpars[ipars] * gc->gcaltpl[ipars*nq+iq]);
    }
    ipars++;
    for(ig=1;ig<gc->ngcal[im1];ig++) {
      for(iq=0;iq<nq;iq++) {
        gc->tpl[im1*nq + iq] *= exp(locpars[ipars] * gc->gcaltpl[ipars*nq+iq]);
       
      }
      ipars++;
    }
  }

  comp_gcal_apply_update(m,SC->nq,rq,gc->tpl,gc->nell,gc->tpll,gc->bins);
}

void comp_gcal_free(void** data) {
  SmicaComp *SC;
  SC_gcal *gc;
  
  SC = *data;
  gc = SC->data;
  free(gc->ngcal);
  free(gc->gcaltpl);
  free(gc->tpl);
  if (gc->nell>0) {
    free(gc->tpll);
    free(gc->bins);
  }
  free(gc);

  free(SC);
  *data = NULL;
}

SmicaComp* comp_gcal2_init(int q,int m, int npar, int *im, int *jm, double* tpl, error **err) {
 SC_gcal2 *gc;
  SmicaComp *SC;

  gc = malloc_err(sizeof(SC_gcal2),err);
  forwardError(*err,__LINE__,NULL);

  gc->npar = npar;

  gc->rtpl = malloc_err(sizeof(double)*npar*q,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->rtpl,tpl,sizeof(double)*npar*q);

  gc->rqbuf = malloc_err(sizeof(double)*m*m*q,err);
  forwardError(*err,__LINE__,NULL);
  
  gc->im = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  gc->jm = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  
  SC = alloc_SC(npar,q,m,gc, &comp_gcal2_update, &comp_gcal2_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;

}


void comp_gcal2_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_gcal2 *gc;
  int ipar;
  int m2, q,m;
  int im,jm,iq,z;
  double *tpl;

  SC = data;
  gc = SC->data;

  m = SC->m;
  m2 = SC->m*SC->m;
  m = SC->m;
  q = SC->nq;
  
  memset(gc->rqbuf,0,sizeof(double)*m2*q);

  for(ipar=0;ipar<gc->npar;ipar++) {
    im = gc->im[ipar];
    jm = gc->jm[ipar];
    tpl = &(gc->rtpl[ipar*q]);
    for(iq=0;iq<q;iq++) {
      gc->rqbuf[iq*m2+im*m+jm] += 2*locpars[ipar] * tpl[iq];  
      gc->rqbuf[iq*m2+jm*m+im] = gc->rqbuf[iq*m2+im*m+jm];
    }    
  }
  for(z=0;z<m2*q;z++) {
    rq[z] *= exp(gc->rqbuf[z]);
  }
}

void comp_gcal2_free(void** data) {
  SmicaComp *SC;
  SC_gcal2 *gc;
  
  SC = *data;
  gc = SC->data;
  free(gc->im);
  free(gc->jm);
  free(gc->rtpl);
  free(gc->rqbuf);
  free(gc);

  free(SC);
  *data = NULL;
}


SmicaComp* comp_calTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,double *w,int *other, error **err ) {
  SC_calTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  gc = malloc_err(sizeof(SC_calTP),err);
  forwardError(*err,__LINE__,NULL);

  gc->npar = npar;


  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  gc->mT = mT;
  gc->mP = mP;
  gc->TEB[0] = TEB[0];
  gc->TEB[1] = TEB[1];
  gc->TEB[2] = TEB[2];

  gc->calvec = malloc_err(sizeof(double)*m,err);
  forwardError(*err,__LINE__,NULL);
  
  for (i=0;i<m;i++) {
    gc->calvec[i]=1;
  }

  
  gc->im = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->im,im,sizeof(int)*npar);

  gc->w = malloc_err(sizeof(double)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->w,w,sizeof(double)*m*m*2);

  gc->other = malloc_err(sizeof(int)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->other,other,sizeof(int)*m*m*2);

  SC = alloc_SC(npar,q,m,gc, &comp_calTP_update, &comp_calTP_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;

}

void comp_calTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_calTP *gc;
  int i,iq,im1,im2,m,m2;

  SC = data;
  gc = SC->data;

  for(i=0;i<SC->ndim;i++) {
    int im;
    im = gc->im[i];
    gc->calvec[im] = exp(locpars[i]);
    if (gc->TEB[2]==1 && im>gc->mT) {
      gc->calvec[im+gc->mP*gc->TEB[1]] = gc->calvec[im];
    }
  }
  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    int iqo;
    iqo = iq*m2;
    for(im1=0;im1<SC->m;im1++) {
      int imo;
      double k1;
      k1 = gc->calvec[im1];
      imo = iqo+im1*m;
      for(im2=im1;im2<SC->m;im2++) {
        int mpos,im1_prime,im2_prime;
        double w,w_prime;
        mpos = (im1*m+im2)*2;
        //_DEBUGHERE_("%d %d %d %g %g %g",iq,im1,im2,k1,gc->calvec[im2],rq[imo+im2]);
        w = gc->w[mpos];
        w_prime = gc->w[mpos+1];
        im1_prime = gc->other[mpos];
        im2_prime = gc->other[mpos+1];
        //if (iq==0) {
        //  _DEBUGHERE_("%d | %d %d -> %g %g %g %g, %d %d -> %g %g %g %g,",mpos,im1,im2,w,gc->calvec[im1],gc->calvec[im2],w*gc->calvec[im1]*gc->calvec[im2],im1_prime,im2_prime,w_prime,gc->calvec[im1_prime],gc->calvec[im2_prime],w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime]);  
        //  _DEBUGHERE_("%g %g",w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime],gc->calvec[im1]*gc->calvec[im2]);
        // }
        rq[imo+im2] *= w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime];
        //rq[imo+im2] *= gc->calvec[im1]*gc->calvec[im2];
      } 
    }
  }
}

void comp_calTP_free(void** data) {
  SmicaComp *SC;
  SC_calTP *gc;
  
  SC = *data;
  gc = SC->data;

  
  free(gc->im);
  free(gc->calvec);
  free(gc->w);
  free(gc->other);
  free(gc);

  free(SC);
  
  *data = NULL;
}


void comp_beamTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_beamTP *gc;
  int t,iq,im1,im2,m,m2,neigen,offm,offq;
  double cal;

  SC = data;
  gc = SC->data;
  memcpy(gc->pars+1,locpars,sizeof(double)*SC->ndim);

  m = SC->m;
  m2 = m*m;
  neigen = gc->neigen;
  
  for(iq=0;iq<SC->nq;iq++) {
    for(im1=0;im1<m;im1++) {
      for(im2=0;im2<m;im2++) {
        cal = 0;
        offm = im1*m+im2;
        offq = offm*m;
        for(t=0;t<neigen;t++) {
          cal += gc->pars[gc->im[offm*neigen + t]] * gc->modes[offq*neigen + t];
        }
        rq[offq] *= exp(cal);
      }
    }
  }
}

void comp_beamTP_free(void** data) {
  SmicaComp *SC;
  SC_beamTP *gc;
  
  SC = *data;
  gc = SC->data;

  
  free(gc->im);
  free(gc->pars);
  free(gc->modes);
  free(gc);

  free(SC);
  
  *data = NULL;
}

SmicaComp* comp_beamTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,int neigen, double *modes,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  gc = malloc_err(sizeof(SC_beamTP),err);
  forwardError(*err,__LINE__,NULL);

  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  gc->pars = malloc_err(sizeof(double)*(npar+1),err);
  forwardError(*err,__LINE__,NULL);
  gc->pars[0] = 0;

  
  gc->neigen = neigen;
  gc->modes = malloc_err(sizeof(double)*neigen*m*m*q,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->modes,modes,sizeof(double)*neigen*m*m*q);  
  
  gc->im = malloc_err(sizeof(int)*neigen*m*m,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->im,im,sizeof(int)*neigen*m*m);

  SC = alloc_SC(npar,q,m,gc, &comp_beamTP_update, &comp_beamTP_free,err);
  forwardError(*err,__LINE__,NULL);

  return SC;

}
