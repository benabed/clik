from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol_template
from clik.parametric import norename,rename_machine,rename_replace

cdef extern c_parametric *ffp8_gcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ffp8_gal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)



cdef class ffp8_gcib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ffp8_gcib_init;
    self.template_name = "cib_ffp8_100_353.dat"
    self.plugin_name = "ffp8"    

cib_ffp8 = rename_machine(ffp8_gcib,{},rename_replace("gib","cib"))
cib_ffp8_v2 = rename_machine(ffp8_gcib,{},rename_replace("gib","cib"),data_file="cib_ffp8_v2_100_353.dat")

cdef class ffp8_gal(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> ffp8_gal_init;
    self.template_name = "dust_ffp8_100_353_TEB.dat"
    self.plugin_name = "ffp8"


def rename_ffp8_dust_TT(v,rups):
  renint={ "gal545_A_100":"A_ffp8_gal_100_100_TT",
           "gal545_A_143":"A_ffp8_gal_143_143_TT",
           "gal545_A_217":"A_ffp8_gal_217_217_TT",
           "gal545_A_143_217":"A_ffp8_gal_143_217_TT"}
  if v in renint:
    rups[v] = renint[v]



dust_ffp8_TT = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"200"},rename_ffp8_dust_TT)
dust_ffp8_v2_TT = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"200"},rename_ffp8_dust_TT,data_file="dust_ffp8_v2_100_353_TEB.dat")


def rename_ffp8_dust_TE(v,rups):
  renint={ "galf_TE_A_100":"A_ffp8_gal_100_100_TE",
           "galf_TE_A_100_143":"A_ffp8_gal_100_143_TE",
           "galf_TE_A_100_217":"A_ffp8_gal_100_217_TE",
           "galf_TE_A_143":"A_ffp8_gal_143_143_TE",
           "galf_TE_A_217":"A_ffp8_gal_217_217_TE",
           "galf_TE_A_143_217":"A_ffp8_gal_143_217_TE"}
  if v in renint:
    rups[v] = renint[v]

def rename_ffp8_dust_EE(v,rups):
  renint={ "galf_EE_A_100":"A_ffp8_gal_100_100_EE",
           "galf_EE_A_100_143":"A_ffp8_gal_100_143_EE",
           "galf_EE_A_100_217":"A_ffp8_gal_100_217_EE",
           "galf_EE_A_143":"A_ffp8_gal_143_143_EE",
           "galf_EE_A_217":"A_ffp8_gal_217_217_EE",
           "galf_EE_A_143_217":"A_ffp8_gal_143_217_EE"}
  if v in renint:
    rups[v] = renint[v]


dust_ffp8_TE = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"500"},rename_ffp8_dust_TE)
dust_ffp8_EE = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"500"},rename_ffp8_dust_EE)

dust_ffp8_v2_TE = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"500"},rename_ffp8_dust_TE,data_file="dust_ffp8_v2_100_353_TEB.dat")
dust_ffp8_v2_EE = rename_machine(ffp8_gal,{"ffp8_gal_l_pivot":"500"},rename_ffp8_dust_EE,data_file="dust_ffp8_v2_100_353_TEB.dat")
 
      
component_list = ["cib_ffp8","dust_ffp8_TT","dust_ffp8_EE","dust_ffp8_TE","cib_ffp8_v2","dust_ffp8_v2_TT","dust_ffp8_v2_EE","dust_ffp8_v2_TE","ffp8_gal","ffp8_gcib"]
