#try to support many flavours of lapack
import autoinstall_lib as atl
from waflib import Logs
import os.path as osp


def options(ctx):
  atl.add_lib_option("lapack",ctx,install=False)
  ctx.add_option("--lapack_mkl",action="store",default="",help="if lapack is mkl, location of the mkl install")

def configure(ctx):
  # no veclib !
  if "vecLib" in ctx.options.lapack_lib or "veclib" in ctx.options.lapack_include:
    #must be darwin !
    if ctx.env.has_lapack:
      raise Exception("Accelerate framework is not supported at this time")
    else:
      Logs.pprint("RED","Apple vecLib not supported")
      return
  
  lapack_extradefs = ["HAS_LAPACK"]
  lapack_libs = ["BLAS","LAPACK"]
  lapack_includes = ["lapack.h","blas.h"]

  if "mkl" in ctx.options.lapack_lib.lower() or "mkl" in ctx.options.lapack_include.lower() or "mkl" in ctx.options.lapack_link or ctx.options.lapack_mkl:
    ctx.env.mkl = True
    lapack_extradefs += ["HAS_MKL"]
    lapack_includes = ["mkl_lapack.h","mkl_blas.h"]
    if ctx.options.lapack_mkl:
      if "framework" in ctx.options.lapack_mkl.lower():
        # guess we are on macosx
        # get the path of the framework
        if ctx.options.lapack_mkl[-1] == "/":
          fpath,fname = osp.split(ctx.options.lapack_mkl[:-1])
        else:
          fpath,fname = osp.split(ctx.options.lapack_mkl)
        fname = fname.split(".")[0]
        ctx.options.lapack_include =  ctx.options.lapack_mkl+"/Headers"
        ctx.options.lapack_lib =  ctx.options.lapack_mkl+"/Libraries/universal"
        if ctx.options.lapack_link=="":
          ctx.options.lapack_link = "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
      else:
        # assume it's 10 on linux
        # check whether it's 10.3
        libsuffix="/lib/em64t"
        if ctx.options.lapack_link=="":
          ctx.options.lapack_link = "-lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread -lmkl_def"
        if osp.exists(ctx.options.lapack_mkl+"/lib/intel64"):
          libsuffix="/lib/intel64"
          ctx.options.lapack_link = "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread"
        ctx.options.lapack_include=ctx.options.lapack_mkl+"/include"
        ctx.options.lapack_lib=ctx.options.lapack_mkl+libsuffix
        
  atl.conf_lib(ctx,"lapack",lapack_libs,"dpotrf",lapack_includes,defines=lapack_extradefs)
