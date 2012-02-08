cimport numpy as nm
import numpy as nm
nm.import_array()
cimport stdlib as stdlib
cimport stdio as stdio

cdef extern from "errorlist.h":
  ctypedef struct error:
    pass
    
  void stringError(char* str, error *err)
  int getErrorValue(error* err)
  int isError(error *err) 
  void purgeError(error **err) 
  void printError(void* flog,error* err)

class CError(Exception):
  def __init__(self,val,str):
    self.val=val
    self.comment=str
  def __str__(self):
    return self.comment.strip()

cdef doError(error **err):
  cdef char estr[10000]
  if (isError(err[0])):
    stringError(estr,err[0])
    er=CError(getErrorValue(err[0]),estr)
    purgeError(err)
    
    return er
  return None


cdef extern from "clik_parametric.h":
  ctypedef struct c_parametric "parametric":
    int lmin,lmax,ndet,nfreq,nvar

  c_parametric *parametric_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  void parametric_free(void** pegl)
  void parametric_compute(c_parametric *egl, double *pars, double* Rq, double *dRq, error **err)

  c_parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

  c_parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *radiogal_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)


cdef class parametric:
  cdef c_parametric* celf
  cdef int nell

  def __init__(self,detlist,vars,lmin,lmax,defs={}):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef error *_err,**err
    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)
    
    for i in range(ndet):
      p_detlist[i] = detlist[i]

    i = 0
    for k,v in defs.items():
      defkey[i] = k
      defvalue[i] = v
      i+=1
    
    nvar = len(vars)
    for i in nvar:
      key[i] = vars[i]

    self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    self.nell = lmax+1-lmin

  def __call__(self,pars,derivatives=False):
    cdef error *_err,**err
    cdef double *_drq,*_rq
      
    
    if len(pars)!=self.celf.nvar:
      raise Exception("Bad shape (expecting (%d) got (%d))"%(self.celf.nvar,len(pars)))
    rq = nm.zeros((self.nell,self.celf.ndet,self.celf.ndet),dtype=nm.double)
    _rq = <double*> nm.PyArray_DATA(rq)
    if derivatives:
      drq = nm.zeros((self.celf.nvar,self.nell,self.celf.ndet,self.celf.ndet),dtype=nm.double)
      _drq = <double*> nm.PyArray_DATA(drq)
      
    else:
      _drq = NULL
    pars_proxy=nm.PyArray_ContiguousFromAny(pars,nm.NPY_DOUBLE,1,1)
    _err = NULL
    err = &_err
    parametric_compute(self.celf,  <double*> nm.PyArray_DATA(pars_proxy), _rq,_drq, err)
    er=doError(err)
    if er:
      raise er
    if derivatives:
      return rq,drq
    return rq
    
  def __dealloc__(self):
    if self.celf!=NULL:
      parametric_free(<void**>&(self.celf))

cdef class powerlaw(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={}):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef error *_err,**err
    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)
    
    for i in range(ndet):
      p_detlist[i] = detlist[i]


    i = 0
    for k,v in defs.items():
      defkey[i] = k
      defvalue[i] = v
      i+=1
    
    nvar = len(vars)
    for i in range(nvar):
      key[i] = vars[i]

    self.celf = powerlaw_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    self.nell = lmax+1-lmin

cdef class powerlaw_free_emissivity(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={}):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef error *_err,**err
    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)
    
    for i in range(ndet):
      p_detlist[i] = detlist[i]

    i = 0
    for k,v in defs.items():
      defkey[i] = k
      defvalue[i] = v
      i+=1
    
    nvar = len(vars)
    for i in range(nvar):
      key[i] = vars[i]

    self.celf = powerlaw_free_emissivity_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    self.nell = lmax+1-lmin

cdef class radiogal(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={}):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef error *_err,**err
    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)
    
    for i in range(ndet):
      p_detlist[i] = detlist[i]

    i = 0
    for k,v in defs.items():
      defkey[i] = k
      defvalue[i] = v
      i+=1
    
    nvar = len(vars)
    for i in range(nvar):
      key[i] = vars[i]

    self.celf = radiogal_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    self.nell = lmax+1-lmin

