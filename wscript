from waflib import Logs
import sys
import os.path as osp

from waflib import Logs
import sys
import os.path as osp

sys.path+=["waf_tools"]
import autoinstall_lib as atl

def options(ctx):
  ctx.load("local_install","waf_tools")
  ctx.load("try_icc","waf_tools")
  ctx.load("try_ifort","waf_tools")
  ctx.load("mbits","waf_tools")
  ctx.load("autoinstall_gsl","waf_tools")
  ctx.load("autoinstall_hdf5","waf_tools")
  ctx.load("any_lapack","waf_tools")
  ctx.load("pmclib","waf_tools")
  ctx.load("chealpix","waf_tools")
  ctx.load("pmclib","waf_tools")
  ctx.load("python")
  try:
    import waflib.Configure
    if not osp.exists(osp.join(osp.split(waflib.__file__)[0],"extras/cython.py")):
      waflib.Configure.download_tool("cython",ctx=ctx)  
    ctx.load("cython",dowload=True)
  except Exception,e:
    pass
  options_h5py(ctx)
  ctx.add_option("--no_pytools",action="store_true",default=False,help="do not build the python tools")
  ctx.add_option("--no_bopix",action="store_true",default=False,help="do not build the python tools")

  
def configure(ctx):
  
  ctx.load("try_icc","waf_tools")
  ctx.load("mbits","waf_tools")
  ctx.load("osx_shlib","waf_tools")
  ctx.load("try_ifort","waf_tools")
  
  ctx.load("local_install","waf_tools")
  
  # dl
  ctx.check_cc(lib="dl",mandatory=1,uselib_store="dl")
  ctx.check_cc(lib="dl",mandatory=0,defines=["HAS_RTLD_DEFAULT"],fragment="#include <dlfcn.h> \nint main() {void* tt = RTLD_DEFAULT;}",msg="checking for RTLD_DEFAULT in dl",uselib_store="dl")
  
  # rpath
  ctx.env.append_value("RPATH",ctx.env.PREFIX+"/lib")

  #configure pmc
  ctx.env.has_pmc = False
  ctx.load("pmclib","waf_tools")

  if (not ctx.env.has_pmc) or ("hdf5" not in ctx.env.LIB_pmc):
    #configure hdf5
    ctx.env.has_hdf5 = True
    ctx.load("autoinstall_hdf5","waf_tools")
  
  if (not ctx.env.has_pmc) or ("HAS_LAPACK" not in ctx.env.DEFINES_pmc):
    #configure lapack
    ctx.env.has_lapack = True
    ctx.load("any_lapack","waf_tools")

  if (not ctx.env.has_pmc) or ("gsl" not in ctx.env.LIB_pmc):
    # configure gsl
    ctx.env.has_gsl = True
    ctx.load("autoinstall_gsl","waf_tools")

  #configure chealpix
  ctx.env.has_chealpix = True
  ctx.load("chealpix","waf_tools")
  
  

  if not ctx.options.no_pytools:
    ctx.load("python")
    if ctx.env.PYTHON[0]!=sys.executable:
      from waflib.Logs import warn
      warn("reverting to current executable")
      ctx.env.PYTHON[0]=sys.executable
      import os
      import os.path as osp
      os.environ["PATH"]=":".join(set(os.environ["PATH"].split(":")+[osp.dirname(sys.executable)]))
    try:
      ctx.check_python_headers()
      # remove unwanted flags for darwin
      _remove_arch(ctx,"CFLAGS_PYEXT")
      _remove_arch(ctx,"LINKFLAGS_PYEXT")
      import numpy
      ctx.env.append_value("INCLUDES_PYEXT",numpy.get_include())
      from distutils.sysconfig import get_config_var
      ctx.env.LIB_PYEXT=ctx.env.LIB_PYEMBED     
      ctx.env.LIBPATH_PYEXT=ctx.env.LIBPATH_PYEMBED+[get_config_var("LIBPL")]     
      import waflib.Configure
      #waflib.Configure.download_tool("cython",ctx=ctx)
      ctx.check_tool("cython")
      configure_h5py(ctx)
    except Exception,e:
      ctx.options.no_pytools = True
      Logs.pprint("BLUE","No suitable python distribution found")
      Logs.pprint("BLUE","Cause : '%s'"%e)
      Logs.pprint("BLUE","Compilation will continue without it (but I strongly advise that you install it)")
      
      
  print "configure ok\n\nrun './waf build install' now !"
  
  
def build(ctx):
  ctx.recurse("src")
  if not ctx.options.no_pytools:
    ctx.recurse("python")
  ctx.add_post_fun(post)

def _remove_arch(ctx,evn):
  if sys.platform.lower()=="darwin":
    cflags_pyext = getattr(ctx.env,evn)
    cflags_pyext_new = []
    inarch = 0
    for cf in cflags_pyext:
      #print cf,inarch
      if inarch == 1:
        inarch = 0
        continue
      if cf == "-arch":
        inarch = 1
        continue
      cflags_pyext_new +=[cf]
    setattr(ctx.env,evn,cflags_pyext_new)
def _prepare_src(ctx):
  import shutil
  import os
  
  # deals with .h
  [shutil.copy("../include/target/%s"%ff,"src") for ff in [
    "aplowly.h",
    "erfinv.h",
    "fowly.h",
    "lklbs.h",
    "lowly_common.h",
    "lowly_pol.h",
    "lowly.h",
    "smica.h"
    ]]
  
  # deals with .c
  [shutil.copy("../src/target/%s"%ff,"src") for ff in [
    "aplowly.c",
    "erfinv.c",
    "fowly.c",
    "lklbs.c",
    "lowly_common.c",
    "lowly_pol.c",
    "lowly.c",
    "smica.c",
    ]]


def dist(ctx):
  import re
  try:
    _prepare_src(ctx)
  except Exception,e:
    pass
  ctx.base_name = 'clik-1.5'
  res = ctx.cmd_and_log("cd ..;svn log -r BASE")
  svnversion = re.findall("(r\d+)",res)[0]
  f=open("svnversion","w")
  print >>f,svnversion
  f.close()
  ctx.files = ctx.path.ant_glob("svnversion waf wscript examples/*.par examples/*.dat **/wscript python/**/*.py python/**/*.pyx src/* src/minipmc/* waf_tools/*.py" )
  
def post(ctx):
  import shutil
  from waflib import Utils
  import os
  if ctx.cmd == 'install':
    # install the module file. This is a cheap trick... grml
    shutil.copy('build/clik.mod',ctx.env.PREFIX+'/include/')
    # go around a waf bug which set the wrong chmod to fortran exec
    os.chmod("bin/clik_example_f90",Utils.O755)
def options_h5py(ctx):
  ctx.load("python")
  ctx.add_option("--h5py_install",action="store_true",default="",help="try to install h5py")

def configure_h5py(ctx):
  import waflib.Logs
  import os
  from waflib import Errors
  ctx.load("python")
  doit = False
  import sys
  sys.path+=[ctx.env.PYTHONDIR]

  if ctx.options.h5py_install:
    try:
      import h5py
      #raise Exception()
    except Exception,e: 
      doit=True
    if doit:
      import os.path as osp
      import autoinstall_lib as atl
      atl.installsmthg_pre(ctx,"http://h5py.googlecode.com/files/h5py-1.3.1.tar.gz","h5py-1.3.1.tar.gz")
      if ctx.env.INCLUDES_hdf5:
        HDF5_DIR=osp.split(ctx.env.INCLUDES_hdf5[0])[0]
      else:
        fi = ctx.find_file("hdf5.h",ctx.env.INCLUDES_pmc)
        #print fi
        HDF5_DIR=osp.split(osp.split(fi)[0])[0]
      HDF5_API="18"
      #print HDF5_DIR
      try:
        os.makedirs(ctx.env.PYTHONDIR)
      except Exception,e:
        print e
      cmdline =  "cd build/%s; HDF5_DIR=%s HDF5_API=%s PYTHONPATH=%s %s setup.py install --install-lib=%s"%("h5py-1.3.1",HDF5_DIR,HDF5_API,ctx.env.PYTHONDIR,ctx.env.PYTHON[0],ctx.env.PYTHONDIR)
      waflib.Logs.pprint("PINK",cmdline)
      ret = ctx.exec_command(cmdline)
      if ret!=0:
        raise Errors.ConfigurationError("Cannot build h5py")
      eggdir = [v for v in os.listdir(ctx.env.PYTHONDIR) if "h5py" in v][0]
      if eggdir!="h5py":
        mdir = [v for v in os.listdir(osp.join(ctx.env.PYTHONDIR,eggdir)) if "h5py" in v][0]
        import os
        os.symlink(osp.join(ctx.env.PYTHONDIR,eggdir,mdir),osp.join(ctx.env.PYTHONDIR,"h5py"))
  try:
    import h5py
  except Exception,e: 
    if not doit:
      waflib.Logs.pprint("PINK", "You can install automatically h5py using cmdline option --h5py_install")      
    else:
      waflib.Logs.pprint("RED", "Autoinstall h5py has failed !")      
    raise e
