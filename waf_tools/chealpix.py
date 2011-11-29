import autoinstall_lib as atl
from waflib import Logs
from waflib import Utils,Errors
import os.path as osp
    
def options(ctx):
  atl.add_lib_option("healpix",ctx,install=True)
  
def configure(ctx):
  iall = atl.shouldIinstall_all(ctx,"healpix")
  if ctx.options.healpix_install or ctx.options.healpix_forceinstall or iall:
    #print "do install"
    ctx.options.healpix_islocal=True
    ctx.options.healpix_forceinstall=True
  #try:
  #  atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cshlib"])
  #  atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fcshlib"])
  #except:
    #if not ctx.options.healpix_install:
    #  raise
    #else:
    #Logs.pprint("PINK","healpix not found. Try to install it")
    Logs.pprint("PINK","Try to install healpix it")
    ctx.options.healpix_islocal=True
    #atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
    #if not bool(ctx.env.has_cfitsio):
    install_cfitsio(ctx)
    atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
    if not bool(ctx.env.has_cfitsio):
      raise Errors.WafError("Cannot build %s"%"cfitsio")
    install_healpix(ctx)
  atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="",opt_name="healpix",uselib=["cshlib"])
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
  import os
  hpdir = "Healpix_2.20a"
  atl.installsmthg_pre(ctx,"http://sourceforge.net/projects/healpix/files/Healpix_2.20a/Healpix_2.20a_2011Feb09.tar.gz/download","Healpix_2.20a_2011Feb09.tar.gz")
  fpic_c = [vv for vv in ctx.env.CFLAGS_cshlib if "-fpic" in vv.lower()]
  fpic_f90 = [vv for vv in ctx.env.CFLAGS_cshlib if "-fpic" in vv.lower()]
  
  dii={"CC":ctx.env.CC[0],"CFLAGS":" ".join(ctx.env.CCFLAGS+fpic_c),"LIBDIR":ctx.env.LIBDIR,"INCDIR":ctx.env.PREFIX+"/include","FC":ctx.env.FC,"FFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fc_omp+fpic_f90)}
  # if I am here, I found cfitsio
  # could it be somewhere else ?
  #cfitsiopath=""
  #for pth in [ctx.env.LIBDIR,"/usr/local/lib","/usr/lib","/lib"]:
  #  print pth    
  #  print osp.join(pth,ctx.env.cshlib_PATTERN%"cfitsio")
  #  if osp.exists(osp.join(pth,ctx.env.cshlib_PATTERN%"cfitsio")):
  #    cfitsiopath = pth
  #    break
  #if not bool(cfitsiopath):
  #  raise Exception("cannot find cfitsio !")
  #dii["CFITSIOPATH"]=cfitsiopath
  #dii["CFITSIOPATHINC"]=osp.realpath(cfitsiopath+"/../include")
  dii["CFITSIOPATH"]=ctx.env.LIBPATH_cfitsio[0]
  dii["CFITSIOPATHINC"]=ctx.env.INCLUDES_cfitsio[0]
  #print dii

  f=open(osp.join("build",hpdir,"conf.cmd"),"w")
  print >>f,cnf_tmpl%dii
  f.close()
  # prepare a few things
  try:
    os.mkdir(osp.join("build",hpdir,"lib"))
  except Exception,e:
    #print e
    pass
  try:
    os.mkdir(osp.join("build",hpdir,"include"))
  except Exception,e:
    #print e
    pass
  cmdline = "cd build/%s; ./configure <conf.cmd; make"%hpdir
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%"Healpix")
  import shutil,os
  for fi in os.listdir(osp.join("build",hpdir,"lib")):
    #print "copy",osp.join("build",hpdir,"lib",fi),osp.join(ctx.env.LIBDIR,fi)
    shutil.copyfile(osp.join("build",hpdir,"lib",fi),osp.join(ctx.env.LIBDIR,fi))
  for fi in os.listdir(osp.join("build",hpdir,"include")):
    #print "copy",osp.join("build",hpdir,"include",fi),osp.join(ctx.env.PREFIX,"include",fi)
    shutil.copyfile(osp.join("build",hpdir,"include",fi),osp.join(ctx.env.PREFIX,"include",fi))
    
cnf_tmpl="""6
n
2
%(CC)s
-O2 -Wall %(CFLAGS)s


%(CFITSIOPATH)s
%(CFITSIOPATHINC)s
y
3
%(FC)s


-I$(F90_INCDIR) %(FFLAGS)s

%(CC)s
-O %(CFLAGS)s


%(CFITSIOPATH)s


0
"""