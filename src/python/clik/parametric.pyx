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


cdef class parametric:
  
  def __cinit__(self):
    self.initfunc=NULL

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

    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<simple_init>self.initfunc)(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
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
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_init

cdef class powerlaw_free_emissivity(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_init

def get_data_path(plg=None):
  import os
  if "CLIK_DATA" in os.environ:
    res=os.environ["CLIK_DATA"]
    if plg:
      return osp.join(res,plg)
  return ""

cdef class parametric_template(parametric):
  def __cinit__(self):
    self.initfunc = NULL
    self.template_name = ""
    self.plugin_name = ""

  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,data_dir="",data_path="",data_file="",data=None):
    cdef int p_detlist[200]
    cdef char *defkey[200],*defvalue[200],*key[200]
    cdef double *template
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


    if data is None:
      if data_path:
        pth = data_path
      else:
        bpth = get_data_path(self.plugin_name)
        if data_dir:
          bpth = data_dir
        fpth = self.template_name
        if data_file:
          fpth = data_file
        pth = osp.join(bpth,fpth)
      tmp = nm.loadtxt(pth)
    else:
      tmp = nm.array(data)
    template = <double*> nm.PyArray_DATA(tmp)
    
    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<template_init>self.initfunc)(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,template,err)
        
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail)

component_list = []
simple_parametric_list = component_list

def register_plugin(plg):
  import sys
  mlg =__import__("clik."+plg,fromlist=[plg])
  global component_list
  component_list += mlg.component_list
  for cp in mlg.component_list:
    globals()[cp]=getattr(mlg,cp)

import os
plgs = [plg.strip() for plg in os.environ.get("CLIK_PLUGIN","").split(",") if plg.strip()]

for plg in plgs:
  try:
    register_plugin(plg)
  except Exception,e:
    pass

