from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template
from clik.parametric import norename,rename_machine,rename_replace

cdef extern double c_sz_spectrum "sz_spectrum" (double nu, double nu0)
cdef extern c_parametric *cib_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *cibr_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *sz_cib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *sz_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *sz_cib_x_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ncib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *tcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric  *ncibXsz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
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

cdef class ncib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ncib_init;
    self.template_name = "cib_model_100_143_217_353.dat"
    self.plugin_name = "cibsz"

cdef class tcib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> tcib_init;
    self.template_name = "cib_tot_100_353_Jsr-1_GL_2013_10.dat"
    self.plugin_name = "cibsz"

cdef class gcib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> tcib_init;
    self.template_name = "cib_1h_2h_100_353_Jsr-1_GL_2014_2.dat"
    self.plugin_name = "cibsz"

cdef class ncibXsz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ncibXsz_init;
    self.template_name = ["ratio_cib_model_100_143_217_353.dat","sz_x_cib_template.dat[1]"]
    self.plugin_name = "cibsz"

cdef class sz_cib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_cib_init;
    self.template_name = ["tsz_143_eps0.50.dat[1]","sz_x_cib_template.dat[1]"]
    self.plugin_name = "cibsz"

cdef class sz_x(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_x_init;
    self.template_name = ["tsz_143_eps0.50.dat[1]","sz_x_cib_template.dat[1]"]
    self.plugin_name = "cibsz"

cdef class sz_cib_x(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_cib_x_init;
    self.template_name = ["tsz_143_eps0.50.dat[1]","sz_x_cib_template.dat[1]"]
    self.plugin_name = "cibsz"

cdef class cib_x(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cib_x_init;
    

cdef class cib(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cib_init;

cdef class cibr(parametric):
  def __cinit__(self):
    self.initfunc = <void*> cibr_init;

    

cib_1h_2h = rename_machine(gcib,{},rename_replace("gcib","cib"),data_file="cib_1h_2h_100_353_Jsr-1_GL_2014_2.dat")
cib_1h_2h_SN = rename_machine(gcib,{},rename_replace("gcib","cib"),data_file="cib_tot_100_353_Jsr-1_GL_2013_10.dat")
cibsz_1h_2h = rename_machine(gcib,{},rename_replace("gcib","cib"),data_file="cibsz_1h_2h_100_353_Jsr-1_GL_2014_2.dat")



      
component_list = ["cibsz_1h_2h","cib_1h_2h","cib_1h_2h_SN","gcib","ncib","cib","cibr","sz","sz_cib","sz_x","cib_x","sz_cib_x","ksz","ncibXsz","tcib"]
