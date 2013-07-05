from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template

cdef extern c_parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pw_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_TE_init

component_list = ["pw_TE"]