import autoinstall_lib as atl
from waflib import Logs
import os.path as osp
    
def options(ctx):
  atl.add_lib_option("healpix",ctx)
  
def configure(ctx):
  atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix")
  atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90")

