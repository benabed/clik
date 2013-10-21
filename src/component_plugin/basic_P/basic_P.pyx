from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template

cdef extern c_parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *pw_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pw_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_TE_init

cdef class pw_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_EE_init


cdef class gal_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_TE_init

cdef class gal_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_EE_init

component_list = ["pw_TE","pw_EE","gal_EE","gal_TE"]
