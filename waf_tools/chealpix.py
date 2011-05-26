import autoinstall_lib as atl
from waflib import Logs
from waflib import Utils,Errors
import os.path as osp
    
def options(ctx):
  atl.add_lib_option("healpix",ctx,install=True)
  
def configure(ctx):
  if ctx.options.healpix_install:
    ctx.options.healpix_islocal=True
  try:
    atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cshlib"])
    atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fcshlib"])
  except:
    if not ctx.options.healpix_install:
      raise
    else:
      Logs.pprint("PINK","healpix not found. Try to install it")
      ctx.options.healpix_islocal=True
      atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
      if not bool(ctx.env.has_cfitsio):
        install_cfitsio(ctx)
        atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
        if not bool(ctx.env.has_cfitsio):
          raise Errors.WafError("Cannot build %s"%"cfitsio")
      install_healpix(ctx)
      atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cshlib"])
      atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fcshlib"])
      
def install_cfitsio(ctx):
  atl.installsmthg_pre(ctx,"ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3280.tar.gz","cfitsio3280.tar.gz")
  CCMACRO = "\"gcc %s\""%ctx.env.mopt
  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
  CPPMACRO = "CPP=\"gcc -E\" CXXCPP=\"g++ -E\" "
  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make;make shared;make install"%("cfitsio",ctx.env.mprefix,"",CCMACRO, CPPMACRO)
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%"cfitsio")
    
def install_healpix(ctx):
  hpdir = "Healpix_2.20a"
  atl.installsmthg_pre(ctx,"http://sourceforge.net/projects/healpix/files/Healpix_2.20a/Healpix_2.20a_2011Feb09.tar.gz/download","Healpix_2.20a_2011Feb09.tar.gz")
  dii={"CC":ctx.env.CC[0],"CFLAGS":" ".join(ctx.env.CCFLAGS+ctx.env.CFLAGS_cshlib),"LIBDIR":ctx.env.LIBDIR,"INCDIR":ctx.env.PREFIX+"/include","FC":ctx.env.FC,"FFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fc_omp+ctx.env.FCFLAGS_fcshlib)}
  f=open(osp.join("build",hpdir,"conf.cmd"),"w")
  print >>f,cnf_tmpl%dii
  f.close()
  cmdline = "cd build/%s; ./configure <conf.cmd; make"%hpdir
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%"Healpix")
  import shutil,os
  for fi in os.listdir(osp.join("build",hpdir,"lib")):
    shutil.copyfile(osp.join("build",hpdir,"lib",fi),osp.join(ctx.env.LIBDIR,fi))
  for fi in os.listdir(osp.join("build",hpdir,"include")):
    shutil.copyfile(osp.join("build",hpdir,"include",fi),osp.join(ctx.env.PREFIX,"include",fi))
    
cnf_tmpl="""2
y
%(CC)s
-O2 -Wall %(CFLAGS)s


%(LIBDIR)s
%(INCDIR)s
y
n
3
%(FC)s


-I$(F90_INCDIR) %(FFLAGS)s

%(CC)s
-O %(CFLAGS)s


%(LIBDIR)s


0
"""