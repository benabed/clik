/*
 *  clik_fortran.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "clik.h"
#ifndef _clik_FORTRAN_
#define _clik_FORTRAN_


// in each of the following functions, if err is set to NULL the code exit as soon as it encounters an error (with a call to exit, bad for you if you're using mpi...)

// initialize the planck likelihood from an hdf file
#ifdef ADD0US
void fortran_clik_init(long* pself,char* hdffilepath,int* fpathlen) {
#elseif ADD2US
void fortran_clik_init__(long* pself,char* hdffilepath,int* fpathlen) {
#else
void fortran_clik_init_(long* pself,char* hdffilepath,int* fpathlen) {
#endif
  clik_object* self;
  hdffilepath[*fpathlen]='\0';
  self = clik_init(hdffilepath,NULL);
  *pself = self; 

}

// retrieve the list of cls as a list of flags 
// order is  TT EE BB TE TB EB
// for example for a likelihood acting on TT 
// has_cl = 1 0 0 0 0 0
#ifdef ADD0US
void fortran_clik_get_has_cl(long* pself, int* has_cl) {
#elseif ADD2US
void fortran_clik_get_has_cl__(long* pself, int* has_cl) {
#else
void fortran_clik_get_has_cl_(long* pself, int* has_cl) {
#endif
  clik_object* self;
  self = *pself;
  clik_get_has_cl(self,has_cl,NULL);

}

// retrieve the number of extra parameters (needed to allocate
// on fortran side)

#ifdef ADD0US
void fortran_clik_get_extra_parameter_number(long* pself, int* numnames) {
#elseif ADD2US
void fortran_clik_get_extra_parameter_number__(long* pself, int* numnames) {
#else
void fortran_clik_get_extra_parameter_number_(long* pself, int* numnames) {
#endif

  clik_object* self;
  parname* names;
  self = *pself;
  *numnames = clik_get_extra_parameter_names(self,&names,NULL);

}

// retrieve the names of extra parameters
#ifdef ADD0US
void fortran_clik_get_extra_parameter_names(long* pself, char* names) {
#elseif ADD2US
void fortran_clik_get_extra_parameter_names__(long* pself, char* names) {
#else
void fortran_clik_get_extra_parameter_names_(long* pself, char* names) {
#endif  
  clik_object* self;
  int i,ii;
  int numnames;
  parname *pnames;
  self = *pself;
  numnames = clik_get_extra_parameter_names(self,&pnames,NULL);

  // Copy parameter names in fortran character array
  for (i=0;i<numnames;i++) {
    memset(&names[i*_pn_size],' ',sizeof(char)*256);
    sprintf(&names[i*_pn_size],"%s",pnames[i]);
  }
  // Get rid of pnames
  free(pnames);

}

// retrieve the lmax for each power spectrum
// -1 --> no cl
// order is TT EE BB TE TB EB
// for example for a likelihood acting ont TT EE TE with same lmax 2000 
// lmax = 2000 2000 -1 2000 -1 -1 -1
#ifdef ADD0US
void fortran_get_lmax(long *pself, int* lmax) {
#elseif ADD2US
void fortran_get_lmax__(long *pself, int* lmax) {
#else
void fortran_get_lmax_(long *pself, int* lmax) {
#endif
  clik_object* self;
  self = *pself;
  clik_get_lmax(self,lmax,NULL);

}

// compute a log likelyhood value
// cl_and_pars is order this way
// first the powerspectra from l=0 to l=lmax[cli] (included) in the order 
// TT EE BB TE TB EB. Only the one where lmax[cli]!=-1 have to be included
// then the extra parameters in the order given by clik_get_extra_parameters.
// The power spectra are in microK^2
// for example, for a likelihood acting on TT, EE and TE with 3 extra parameters 
// will expect an array ordered this way
// C_0^TT ... C_lmax[0]^TT C_0^EE ... C_lmax[1]^EE C_0^TE ... C_lmax[3]^T3 extrapar1 extrapar2 extrapar3

#ifdef ADD0US
void fortran_clik_compute(long* pself, double* cl_and_pars, double* lkl) {
#elseif ADD2US
void fortran_clik_compute__(long* pself, double* cl_and_pars, double* lkl) {
#else
void fortran_clik_compute_(long* pself, double* cl_and_pars, double* lkl) {
#endif
  clik_object* self;
  self = *pself;
  *lkl=clik_compute(self,cl_and_pars,NULL);

}

// cleanup
#ifdef ADD0US
void fortran_clik_cleanup(long* pself) {
#elseif ADD2US
void fortran_clik_cleanup__(long* pself) {
#else
void fortran_clik_cleanup_(long* pself) {
#endif
  clik_object* self;
  self = *pself;
  clik_cleanup(&self);
  *pself = self;

}



#endif
