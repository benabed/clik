from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import powerlaw_free_emissivity,rename_machine,rename_replace

cdef extern c_parametric *gal545_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)


cdef class gal_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_TE_init

cdef class gal_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_EE_init

cdef class gal545(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gal545_init;


component_list = ["gal_EE","gal_TE","gal545"]


cdef extern c_parametric *pointsource_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pointsource(parametric):
  def __cinit__(self):
    self.initfunc = <void*>pointsource_init

component_list += ["pointsource"]


cdef extern double c_sz_spectrum "sz_spectrum" (double nu, double nu0)
cdef extern c_parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric  *gibXsz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *gcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

def sz_spectrum(nu,nu0=143.0):
  return c_sz_spectrum(<double>nu,<double>nu0)

cdef class sz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_init;
    self.template_name = "tsz_143_eps0.50.dat[1]"
    self.plugin_name = "cibsz"

cdef class ksz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ksz_init;
    self.template_name = "ksz_fromcamspec.dat"
    self.plugin_name = "cibsz"

cdef class gcib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> gcib_init;
    self.template_name = "cib_1h_2h_100_353_Jsr-1_PS_2014_09.dat"
    self.plugin_name = "cibsz"

cdef class gibXsz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> gibXsz_init;
    self.template_name = ["sz_x_cib_template.dat[1]"]
    self.plugin_name = "cibsz"
    

cib_1h_2h_sept14 = rename_machine(gcib,{},rename_replace("gib","cib"))
sz_color = rename_machine(sz,{"sz_color_143_to_143":"0.975","sz_color_100_to_143":"0.981"},norename)

component_list += ["sz_color","gibXsz","gcib","sz","ksz","cib_1h_2h_sept14"]

cdef extern c_parametric *bleak_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cnoise_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class bleak(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> bleak_init;
    self.template_name = "sky_template_v15_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class cnoise(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> cnoise_init;
    self.template_name = "cnoise_F100_143_217_353_v17.dat"
    self.plugin_name = "systematics"

cnoise_v17 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_v17.dat")
bleak_v15 = rename_machine(bleak,{},norename,data_file="sky_template_v15_F100_143_217_353.dat")    

component_list = ["bleak","cnoise","cnoise_v17","bleak_v15"]
 

