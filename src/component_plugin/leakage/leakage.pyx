from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import norename,rename_machine,rename_replace

cdef extern c_parametric *bleak_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class bleak(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> bleak_init;
    self.template_name = "sky_template_M40404000_F100_143_217_353.dat"
    self.plugin_name = "leakage"

component_list = ["bleak"]
