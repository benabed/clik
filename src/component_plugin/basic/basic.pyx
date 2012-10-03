from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *radiogal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *galactic_component_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern double c_dust_spectrum "dust_spectrum" (double nu, double T_dust, double beta_dust, double nu0)
cdef extern double c_d_dust_spectrum_d_beta_dust "d_dust_spectrum_d_beta_dust" (double nu, double T_dust, double beta_dust, double nu0)
cdef extern double c_d_dust_spectrum_d_T_dust "d_dust_spectrum_d_T_dust" (double nu, double T_dust, double beta_dust, double nu0)
cdef extern double c_non_thermal_spectrum "non_thermal_spectrum" (double nu, double alpha_non_thermal, double nu0)
cdef extern double c_d_non_thermal_spectrum_d_alpha_non_thermal "d_non_thermal_spectrum_d_alpha_non_thermal" (double nu, double alpha_non_thermal, double nu0)
cdef extern double c_dBdT "dBdT" (double nu, double nu0)
cdef extern double c_sz_spectrum "sz_spectrum" (double nu, double nu0)
cdef extern c_parametric *ir_poisson_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *ir_clustered_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *poisson_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *poisson_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_triangle_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tanh_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
#cdef extern c_parametric *constant_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *sz_cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)


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

def sz_spectrum(nu,nu0=143.0):
  return c_sz_spectrum(<double>nu,<double>nu0)

cdef class radiogal(parametric):
  def __cinit__(self):
    self.initfunc = <void*>radiogal_init


cdef class galametric(parametric):
  def __cinit__(self):
    self.initfunc = <void*>galactic_component_init


cdef class ir_poisson(parametric):
  def __cinit__(self):
    self.initfunc = <void*>ir_poisson_init

cdef class ir_clustered(parametric):
  def __cinit__(self):
    self.initfunc = <void*>ir_clustered_init

cdef class poisson_tensor_bydet(parametric):
  def __cinit__(self):
    self.initfunc = <void*>poisson_tensor_bydet_init


cdef class powerlaw_tensor_bydet(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tensor_bydet_init

cdef class poisson_tensor(parametric):
  def __cinit__(self):
    self.initfunc = <void*>poisson_tensor_init


cdef class powerlaw_tensor(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tensor_init

cdef class powerlaw_triangle(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_triangle_init

cdef class powerlaw_tanh(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tanh_init

#cdef class constant(parametric):
#  def __cinit__(self):
#    self.initfunc = <void*>constant_init
  
cdef class sz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_init;
    self.template_name = "sz.dat"
    self.plugin_name = "basic"
#
#cdef class sz_cib(parametric_template):
#  def __cinit__(self):
#    self.initfunc = <void*> sz_cib_init;
#    self.template_name = "sz_cib.dat"
#    self.plugin_name = "basic"
#
component_list = ["poisson_tensor","powerlaw_tensor","powerlaw_triangle","powerlaw_tanh","poisson_tensor_bydet","powerlaw_tensor_bydet","radiogal","galametric","ir_clustered","ir_poisson","sz","sz_cib"]

