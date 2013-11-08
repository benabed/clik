from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern double c_sz_spectrum "sz_spectrum" (double nu, double nu0)
cdef extern c_parametric *poisson_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *poisson_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_triangle_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tanh_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *pointsource_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *cib_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *cibr_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *sz_cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *sz_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *sz_cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ncib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)



def sz_spectrum(nu,nu0=143.0):
  return c_sz_spectrum(<double>nu,<double>nu0)


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

cdef class pointsource(parametric):
  def __cinit__(self):
    self.initfunc = <void*>pointsource_init

cdef class sz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_init;
    self.template_name = "sz.dat"
    self.plugin_name = "basic"

cdef class ksz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ksz_init;
    self.template_name = "ksz_fromcamspec.dat"
    self.plugin_name = "basic"

cdef class ncib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ncib_init;
    self.template_name = "cib_model_100_143_217_353.dat"
    self.plugin_name = "basic"

cdef class sz_cib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_cib_init;
    self.template_name = "sz_cib.dat"
    self.plugin_name = "basic"

cdef class sz_x(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_x_init;
    self.template_name = "sz_cib.dat"
    self.plugin_name = "basic"

cdef class sz_cib_x(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_cib_x_init;
    self.template_name = "sz_cib.dat"
    self.plugin_name = "basic"

cdef class cib_x(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cib_x_init;
    

cdef class cib(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cib_init;

cdef class cibr(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cibr_init;

    
component_list = ["ncib","cib","cibr","pointsource","poisson_tensor","powerlaw_tensor","powerlaw_triangle","powerlaw_tanh","poisson_tensor_bydet","powerlaw_tensor_bydet","sz","sz_cib","sz_x","cib_x","sz_cib_x","ksz"]

