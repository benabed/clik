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
  options_numpy(ctx)
  options_h5py(ctx)
  options_cython(ctx)
  ctx.add_option("--no_pytools",action="store_true",default=False,help="do not build the python tools")
  ctx.add_option("--no_bopix",action="store_true",default=False,help="do not build the python tools")
  ctx.add_option("--wmap_src",action="store",default="",help="location of wmap likelihood sources")
  ctx.add_option("--wmap_install",action="store_true",default=False,help="download wmap likelihood for me")
  
  
def configure(ctx):
  import os
  import os.path as osp
  ctx.load("try_icc","waf_tools")
  ctx.load("mbits","waf_tools")
  ctx.load("osx_shlib","waf_tools")
  try:
    ctx.load("try_ifort","waf_tools")
    ctx.env.has_f90 = True
  except Exception,e:
    Logs.pprint("BLUE","No fortran compiler found. Will keep on working without it (cause: '%s')"%e)
    ctx.env.has_f90 = False
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
  
  #bopix
  ctx.env.no_bopix = ctx.options.no_bopix or not osp.exists("src/bopix")
  # wmap
  if ctx.options.wmap_install:
    atl.installsmthg_pre(ctx,"http://lambda.gsfc.nasa.gov/data/map/dr4/dcp/wmap_likelihood_sw_v4p1.tar.gz","wmap_likelihood_sw_v4p1.tar.gz","src/")
    ctx.options.wmap_src = "likelihood_v4p1"
  ctx.env.wmap_src =   ctx.options.wmap_src
    
  if not ctx.options.no_pytools:
    ctx.load("python")
    if ctx.env.PYTHON[0]!=sys.executable:
      from waflib.Logs import warn
      warn("reverting to current executable")
      ctx.env.PYTHON[0]=sys.executable
      os.environ["PATH"]=":".join(set(os.environ["PATH"].split(":")+[osp.dirname(sys.executable)]))
    try:
      ctx.check_python_headers()
      # remove unwanted flags for darwin
      _remove_arch(ctx,"CFLAGS_PYEXT")
      _remove_arch(ctx,"LINKFLAGS_PYEXT")
      from distutils.sysconfig import get_config_var
      configure_numpy(ctx)
      configure_h5py(ctx)
      configure_cython(ctx)
      
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
  #ctx.recurse("src/egfs")
  
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
  [os.symlink(osp.realpath("../include/target/%s"%ff),"src/%s"%ff) for ff in [
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
  [os.symlink(osp.realpath("../src/target/%s"%ff),"src/%s"%ff) for ff in [
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
  print "private"
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
  ctx.files = ctx.path.ant_glob("svnversion waf wscript examples/*.par examples/*.dat **/wscript python/**/*.py python/**/*.pyx src/* src/minipmc/* src/bopix/* waf_tools/*.py clik.pdf" )
  
import waflib
class Dist_public(waflib.Scripting.Dist):
  cmd = 'dist_public'
  fun = 'dist_public'
  
def dist_public(ctx):
  print "public"
  import re
  try:
    _prepare_src(ctx)
  except Exception,e:
    pass
  ctx.base_name = 'clik-1.5.public'
  res = ctx.cmd_and_log("cd ..;svn log -r BASE")
  svnversion = re.findall("(r\d+)",res)[0]
  f=open("svnversion","w")
  print >>f,svnversion
  f.close()
  ctx.files = ctx.path.ant_glob("svnversion waf wscript examples/*.par examples/*.dat **/wscript python/**/*.py python/**/*.pyx src/* src/minipmc/* waf_tools/*.py clik.pdf" )
  
def post(ctx):
  import shutil
  from waflib import Utils
  import os
  if ctx.cmd == 'install':
    # install the module file. This is a cheap trick... grml
    shutil.copy('build/clik.mod',ctx.env.PREFIX+'/include/')
    # go around a waf bug which set the wrong chmod to fortran exec
    os.chmod("bin/clik_example_f90",Utils.O755)
    build_env_files(ctx)
    
def build_env_files(ctx):
  import os
  import os.path as ops
  full_libpath = set(ctx.env.LIBPATH_chealpix + ctx.env.LIBPATH_fc_runtime + ctx.env.LIBPATH_gsl + ctx.env.LIBPATH_hdf5 + ctx.env.LIBPATH_healpix_f90 + ctx.env.LIBPATH_lapack)
  
  if osp.basename(os.environ["SHELL"]) in ("csh","tcsh","zsh"):
    name = "clik_profile.csh"
    shebang = "#! /bin/tcsh"
    block_tmpl = """
if !($?%(VAR)s) then
  setenv %(VAR)s %(PATH)s
else
  setenv %(VAR)s %(PATH)s:${%(VAR)s}
endif
"""
  else:
    name = "clik_profile.sh"
    shebang = "#! /bin/sh"
    block_tmpl = """
if [ -z "${%(VAR)s}" ]; then
  %(VAR)s=%(PATH)s
else
  %(VAR)s=%(PATH)s:${%(VAR)s}
fi
export %(VAR)s
"""
  if sys.platform.lower()=="darwin":
    LD_LIB = "DYLD_LIBRARY_PATH"
  else:
    LD_LIB = "LD_LIBRARY_PATH"
  f = open(osp.join(ctx.env.BINDIR,name),"w")
  print >>f,"# this code cannot be run directly"
  print >>f,"# do 'source %s' from your %s shell or put it in your profile"%(osp.join(ctx.env.BINDIR,name),osp.basename(os.environ["SHELL"]))
  print >>f,block_tmpl%{"PATH":ctx.env.BINDIR,"VAR":"PATH"}
  print >>f,block_tmpl%{"PATH":ctx.env.PYTHONDIR,"VAR":"PYTHONPATH"}
  print >>f,block_tmpl%{"PATH":":".join(full_libpath),"VAR":LD_LIB}
  f.close()
  
  print "Use %s to set the environment variables needed by clik"%osp.join(ctx.env.BINDIR,name)
  
def options_h5py(ctx):
  ctx.load("python")
  ctx.add_option("--h5py_install",action="store_true",default="",help="try to install h5py")
def options_numpy(ctx):
  ctx.load("python")
  ctx.add_option("--numpy_install",action="store_true",default="",help="try to install numpy")
def options_cython(ctx):
  ctx.load("python")
  ctx.add_option("--cython_install",action="store_true",default="",help="try to install numpy")

def configure_h5py(ctx):
  import autoinstall_lib as atl
  cmdline=None
  if ctx.options.h5py_install:
    if ctx.env.INCLUDES_hdf5:
      HDF5_DIR=osp.split(ctx.env.INCLUDES_hdf5[0])[0]
    else:
      fi = ctx.find_file("hdf5.h",ctx.env.INCLUDES_pmc)
      #print fi
      HDF5_DIR=osp.split(osp.split(fi)[0])[0]
    HDF5_API="18"
    cmdline =  "cd build/%s; HDF5_DIR=%s HDF5_API=%s PYTHONPATH=%s %s setup.py install --install-lib=%s"%("h5py-1.3.1",HDF5_DIR,HDF5_API,ctx.env.PYTHONDIR,ctx.env.PYTHON[0],ctx.env.PYTHONDIR)

  atl.configure_python_module(ctx,"h5py","http://h5py.googlecode.com/files/h5py-1.3.1.tar.gz","h5py-1.3.1.tar.gz","h5py-1.3.1",cmdline)

def configure_numpy(ctx):
  import autoinstall_lib as atl
  atl.configure_python_module(ctx,"numpy","http://sourceforge.net/projects/numpy/files/NumPy/1.6.0/numpy-1.6.0.tar.gz/download","numpy-1.6.0.tar.gz","numpy-1.6.0")
  import numpy
  ctx.env.append_value("INCLUDES_PYEXT",numpy.get_include())
  
def configure_cython(ctx):
  import autoinstall_lib as atl
  from waflib import Utils
  import os.path as osp
  import os
  
  vv=False
  try:
    # check for cython
    atl.check_python_module(ctx,"cython")
    vv=True
    version_str = "unknown"
    ctx.start_msg("Checking cython version (>0.12)")
    import Cython.Compiler.Version
    version_str = Cython.Compiler.Version.version
    version = [int(v) for v in version_str.split(".")]
    #print version
    assert version[1]>=12
    ctx.end_msg(version_str)
  except Exception,e:
    if vv:
      ctx.end_msg("no (%s)"%version_str,'YELLOW')
    # no cython, install it !
    atl.configure_python_module(ctx,"cython","http://cython.org/release/Cython-0.14.1.tar.gz","Cython-0.14.1.tar.gz","Cython-0.14.1",forceinstall=True)

  try:
    ctx.load("cython")
  except:
    ctx.env.CYTHON=[osp.join(ctx.env.BINDIR,"cython")]
    f=open(osp.join(ctx.env.BINDIR,"cython"))
    cytxt = f.readlines()
    f.close()
    cytxt[1:1] = ["import sys\n","sys.path+=['%s']\n"%str(ctx.env.PYTHONDIR)]
    f=open(osp.join(ctx.env.BINDIR,"cython"),"w")
    f.write("".join(cytxt))
    f.close()
    os.chmod(osp.join(ctx.env.BINDIR,"cython"),Utils.O755)

