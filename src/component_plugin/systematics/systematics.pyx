from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import norename,rename_machine,rename_replace

cdef extern c_parametric *bleak_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cnoise_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *dip_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class bleak(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> bleak_init;
    self.template_name = "sky_template_M40404000_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class cnoise(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> cnoise_init;
    self.template_name = "cnoise_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class dip(parametric):
  def __cinit__(self):
    self.initfunc = <void*> dip_init;

cnoise_gpe = rename_machine(cnoise,{},norename,data_file="cnoise_GPE_F100_143_217_353.dat")
cnoise_t2 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_t2.dat")
cnoise_t3 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_t3.dat")
    
component_list = ["bleak","cnoise","dip","cnoise_gpe","cnoise_t2","cnoise_t3"]
