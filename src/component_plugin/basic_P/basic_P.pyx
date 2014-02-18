from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import powerlaw_free_emissivity,rename_machine

cdef extern c_parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *pw_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef extern c_parametric *powerlaw_free_emissivity_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_free_emissivity_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pw_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_TE_init

cdef class pw_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_EE_init

cdef class powerlaw_free_emissivity_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_EE_init

cdef class powerlaw_free_emissivity_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_TE_init

def galf_P_rename_func(v,rups):
  if v.startswith("galf"):
    rv = v.replace("galf","pwfe")
    rups[v]=rv

galf_TE = rename_machine(powerlaw_free_emissivity_TE,{"galf_TE_l2_norm":"1","galf_TE_l_pivot":"500","ps_TE_index":"-2.4"},galf_P_rename_func)
galf_EE = rename_machine(powerlaw_free_emissivity_EE,{"galf_EE_l2_norm":"1","galf_EE_l_pivot":"500","ps_EE_index":"-2.4"},galf_P_rename_func)

def ps_P_rename_func(v,rups):
  if v.startswith("ps"):
    rv = v.replace("ps","pwfe")
    rups[v]=rv

ps_TE = rename_machine(powerlaw_free_emissivity_TE,{"ps_TE_l2_norm":"1","ps_TE_l_pivot":"3000","ps_TE_index":"0"},ps_P_rename_func)
ps_EE = rename_machine(powerlaw_free_emissivity_EE,{"ps_EE_l2_norm":"1","ps_EE_l_pivot":"3000","ps_EE_index":"0"},ps_P_rename_func)

component_list = ["pw_TE","pw_EE","powerlaw_free_emissivity_EE","powerlaw_free_emissivity_TE","galf_TE","galf_EE","ps_TE","ps_EE"]
