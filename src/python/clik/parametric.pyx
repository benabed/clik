cimport numpy as nm
import numpy as nm
nm.import_array()
cimport stdlib as stdlib
cimport stdio as stdio
import os.path as osp

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
  double parametric_get_default(c_parametric* egl,char *key, error **err)
  double parametric_get_value(c_parametric *egl, char *key, error **err)
  void parametric_dnofail(c_parametric* egl, int vl)


  c_parametric *powerlaw_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *powerlaw_free_emissivity_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *radiogal_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *galactic_component_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  double c_dust_spectrum "dust_spectrum" (double nu, double T_dust, double beta_dust, double nu0)
  double c_d_dust_spectrum_d_beta_dust "d_dust_spectrum_d_beta_dust" (double nu, double T_dust, double beta_dust, double nu0)
  double c_d_dust_spectrum_d_T_dust "d_dust_spectrum_d_T_dust" (double nu, double T_dust, double beta_dust, double nu0)
  double c_non_thermal_spectrum "non_thermal_spectrum" (double nu, double alpha_non_thermal, double nu0)
  double c_d_non_thermal_spectrum_d_alpha_non_thermal "d_non_thermal_spectrum_d_alpha_non_thermal" (double nu, double alpha_non_thermal, double nu0)
  double c_dBdT "dBdT" (double nu, double nu0)
  c_parametric *ir_poisson_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *ir_clustered_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *ir_poisson_pep_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double *rq_poisson_in, error **err)  
  c_parametric *ir_clustered_pep_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double *rq_clustered_in, error **err)

cdef class parametric:
  cdef c_parametric* celf
  cdef int nell
  cdef readonly object varpar,parvalues

  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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
    
    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

  def _post_init(self,detlist,vars,lmin,lmax,defs,dnofail):
    parametric_dnofail(self.celf,int(dnofail))
    prs = vars
    if not dnofail:
      prs = [p for p in vars if self.has_parameter(p)]
      if len(prs)!= len(vars):
        parametric_free(<void**>&(self.celf))
        self.__init__(detlist,prs,lmin,lmax,defs,False)

    
    dv = []
    for p in prs:
      if self.has_parameter(p):
        dv += [self.get_default_value(p)]
      else:
        dv +=[0]

    self.varpar = prs
    self.parvalues = dv
    #self.varpar = OrderedDict(zip(prs,dv))
    self.nell = lmax+1-lmin    

  def get_default_value(self,key):
    cdef error *_err,**err
    _err = NULL
    err = &_err
    
    res = parametric_get_default(self.celf,key, err)
    er=doError(err)
    if er:
      raise er
    return res


  def has_parameter(self,key):
    try:
      self.get_default_value(key)
      return True
    except Exception,e:
      return False

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
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

cdef class powerlaw_free_emissivity(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)


cdef class radiogal(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

cdef class galametric(parametric):
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self.celf = galactic_component_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)


def dust_spectrum(nu,T_dust=18.0,beta_dust=1.8,nu0=143.0):
  return c_dust_spectrum(<double>nu,<double>T_dust,<double>beta_dust,<double>nu0)

def d_dust_spectrum_d_T_dust(nu,T_dust=18.0,beta_dust=1.8,nu0=143.0):
  return c_d_dust_spectrum_d_T_dust(<double>nu,<double>T_dust,<double>beta_dust,<double>nu0)

def d_dust_spectrum_d_beta_dust(nu,T_dust=18.0,beta_dust=1.8,nu0=143.0):
  return c_d_dust_spectrum_d_beta_dust(<double>nu,<double>T_dust,<double>beta_dust,<double>nu0)

def non_thermal_spectrum(nu,alpha_non_thermal=-1.0,nu0=143.0):
  return c_non_thermal_spectrum(<double>nu,<double>alpha_non_thermal,<double>nu0)

def d_non_thermal_spectrum_d_alpha_non_thermal(nu,alpha_non_thermal=-1.0,nu0=143.0):
  return c_d_non_thermal_spectrum_d_alpha_non_thermal(<double>nu,<double>alpha_non_thermal,<double>nu0)

def dBdT(nu,nu0=143.0):
  return c_dBdT(<double>nu,<double>nu0)

cdef class ir_poisson(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self.celf = ir_poisson_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

cdef class ir_clustered(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False):
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

    self.celf = ir_clustered_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

def get_data_path():
  import os
  if "CLIK_DATA" in os.environ:
    return os.environ["CLIK_DATA"]
  return ""

def get_pep_cib_data_path():
  return osp.join(get_data_path(),"pep_cib")

cdef class ir_poisson_pep(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,input_filename=""):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef double *rq_poisson_in
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


    if input_filename:
      tmp = nm.loadtxt(input_filename)
    else:
      tmp = nm.loadtxt(osp.join(get_pep_cib_data_path(),"ir_poisson_pep.dat"))
    rq_poisson_in = <double*> nm.PyArray_DATA(tmp)
        
    self.celf = ir_poisson_pep_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,rq_poisson_in,err)
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

cdef class ir_clustered_pep(parametric):
  
  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,input_filename=""):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef double *rq_clustered_in
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

    if input_filename:
      tmp = nm.loadtxt(input_filename)
    else:
      tmp = nm.loadtxt(osp.join(get_pep_cib_data_path(),"ir_clustered_pep.dat"))
    tmp = nm.loadtxt(input_filename)
    rq_clustered_in = <double*> nm.PyArray_DATA(tmp)
        
    self.celf = ir_clustered_pep_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,rq_clustered_in,err)
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

simple_parametric_list = ["radiogal","ir_clustered","ir_poisson","ir_clustered_pep","ir_poisson_pep","galametric"]
