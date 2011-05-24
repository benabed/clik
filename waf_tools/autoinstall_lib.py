from waflib import Logs
import sys
import os.path as osp

def add_lib_option(libname,opt,default="/usr/local/",install=True):
  opt.add_option("--%s_islocal"%libname,action="store_true",default=False,help="%s has been installed with install%s"%(libname,libname))
  opt.add_option("--%s_prefix"%libname,action="store",default=default,help="%s include/lib path prefix"%libname)
  opt.add_option("--%s_include"%libname,action="store",default="",help="%s include path"%libname)
  opt.add_option("--%s_lib"%libname,action="store",default="",help="%s lib path"%libname)
  opt.add_option("--%s_link"%libname,action="store",default="",help="%s link line"%libname)
  if install:
    opt.add_option("--%s_install"%libname,action="store_true",default="",help="%s try to install"%libname)
    

def noemptylist(li):
  return [ll for ll in li if ll]
  
def opt_to_libpaths(ctx,name):
  if getattr(ctx.options,name+"_islocal",None):
    prefix=ctx.env.localpref
    include=[]
    lib=[]
    link=[]
  else :
    prefix = getattr(ctx.options,name+"_prefix")
    include = getattr(ctx.options,name+"_include").split(":")
    lib = getattr(ctx.options,name+"_lib").split(":")
    link = libsfromlinkline(getattr(ctx.options,name+"_link"))
  include += [osp.join(prefix,"include")]
  lib += [osp.join(prefix,"lib")]
  
  return prefix,include,lib,link
      
def libsfromlinkline(ll):
  return [lb.strip() for lb in ll.split("-l") if lb.strip()]
  
def magic_join(pth,adx):
  if isinstance(pth,str):
    pth = set(pth.split(":"))
  return ":".join([osp.join(pp,adx) for pp in set(pth)])
      
def add_lib(conf,prefix,include,libpath,libname, funcname="",headername="",libs = [], uselib=[],defines=[],frameworkpath=[],framework=[],flagline=""):
  #print conf,prefix,include,libpath,libname, funcname,headername,libs , uselib,defines
  
  for inc in include:
    if inc:
      conf.env.append_value("INCLUDES_%s"%(libname),inc)
  if libs == []:
    libs = [libname]
  if type(uselib)==type(""):
    uselib = [uselib]
  if type(funcname)==type(""):
    funcname=[funcname]
  if type(defines)==type(""):
    defines=[defines]
    
  conf.parse_flags(flagline,uselib=libname)

  conf.check_cc(lib=libs, libpath = noemptylist(libpath),rpath=noemptylist(libpath) ,uselib_store=libname,mandatory=1,uselib=uselib+[libname],defines=defines,frameworkpath=frameworkpath,framework=framework)
  for fnc in funcname:
    conf.check_cc(
      errmsg="failed (check whether lib is compiled in 32 or 64bits)",
      function_name=fnc,header_name=headername,uselib=" ".join([libname]+uselib),mandatory=1,frameworkpath=frameworkpath,framework=framework)

def add_lib_f90(conf,prefix,include,libpath,libname, funcname="",headername="",libs = [], uselib=[],defines=[],frameworkpath=[],framework=[],flagline=""):
  #print conf,prefix,include,libpath,libname, funcname,headername,libs , uselib,defines
  
  for inc in include:
    if inc:
      conf.env.append_value("INCLUDES_%s"%(libname),inc)
  if libs == []:
    libs = [libname]
  if type(uselib)==type(""):
    uselib = [uselib]
  if type(funcname)==type(""):
    funcname=[funcname]
  if type(defines)==type(""):
    defines=[defines]
    
  conf.parse_flags(flagline,uselib=libname)
  # do nothing for now...
  for inc in libpath:
    if inc:
      conf.env.append_value("LIBPATH_%s"%(libname),inc)
      conf.env.append_value("RPATH_%s"%(libname),inc)
  for inc in libs:
    if inc:
      conf.env.append_value("LIB_%s"%(libname),inc)
  return
  #print dir(conf)
  conf.check_fortran(lib=libs, libpath = noemptylist(libpath),rpath=noemptylist(libpath) ,uselib_store=libname,mandatory=1,uselib=uselib+[libname],defines=defines,frameworkpath=frameworkpath,framework=framework)
  for fnc in funcname:
    conf.check_fortran(
      errmsg="failed (check whether lib is compiled in 32 or 64bits)",
      function_name=fnc,header_name=headername,uselib=" ".join([libname]+uselib),mandatory=1,frameworkpath=frameworkpath,framework=framework)



add_lib_dict = {
  "c" : add_lib,
  "f90" : add_lib_f90
}

def conf_lib(ctx,name,_libs,testfunc=[],testinclude=[],add_inc_path=[],defines=[],frameworkpath=[],framework=[],install=False,msg="",uselib=[],flagline="",opt_name="",add_lib_code="c"):
  if not opt_name:
    opt_name=name
  # do install if needed
  if install:
    # first try without install !
    try:
      setattr(ctx.env,"has_"+name,True)
      conf_lib(ctx,name,_libs,testfunc,testinclude,add_inc_path,defines,frameworkpath,framework,False,msg,uselib,flagline,opt_name,add_lib_code)
    except Exception,e:
      Logs.pprint("RED","%s not found, try to install it"%name)
      if getattr(ctx.options,opt_name+"_install"):
        install(ctx)
        setattr(ctx.options,"%s_islocal"%opt_name,1)
  # compute paths
  prefix,include,lib,link = opt_to_libpaths(ctx,opt_name)
  
  # libs to test for
  libs=_libs
  # if a link option is given, it overides the libs
  if link:
    libs = link
    
  # extra includes ?
  if add_inc_path:
    extinc = []
    for adi in add_inc_path:
      extinc += [osp.join(inc,adi) for inc in include]
    include += extinc
  
  try:
    add_lib_dict[add_lib_code](ctx,prefix,include,lib,name,
          testfunc,testinclude,
          libs=libs,uselib=uselib,defines=defines,frameworkpath=frameworkpath,framework=framework,flagline=flagline)

    setattr(ctx.env,"use_%s"%name,name)
    setattr(ctx.env,"has_%s"%name,name)

  except Exception,e:
    ctx.env["INCLUDES_%s"%name]=[]
    if not getattr(ctx.env,"has_"+name,False):
      Logs.pprint("BLUE","Optional %s not found"%name)
      Logs.pprint("BLUE","Compilation will continue without it")
    else:
      Logs.pprint("RED","%s not found"%name)
      Logs.pprint("PINK", "check that %s_prefix or %s_lib and %s_include command line options point toward your %s install"%(name,name,name,name))
      Logs.pprint("PINK", "or check that %s is compiled in %d bit"%(name,{True:64}.get(ctx.options.m64,32)))
      if msg:
        Logs.pprint("PINK", msg)      
      if install:
        Logs.pprint("PINK", "or install automatically using cmdline option --%s_install"%(name))      
      raise e

def installsmthg_pre(ctx,where,what):

  from waflib import Options, Environment,Utils,Errors
  import urllib2
  import re
  import os.path as osp
  import tarfile
  import shutil
  import os
  
  #ctx.env = Environment.Environment(filename="build/c4che/default.cache.py")

  if osp.exists("build/"+what):
    Logs.pprint("PINK","%s already downloaded"%what)
  else:
    Logs.pprint("PINK","download from "+where)
    luaf = urllib2.urlopen(where)
    #if luaf.code!=200 and luaf.code!=None:
    #  raise Utils.WscriptError("Cannot install : %d reported error %d"%(luaf.code,where))
    f=open("build/"+what,"w")
    print >>f,luaf.read(),
    luaf.close()
    f.close()
  tf = tarfile.open("build/"+what)
  #Logs.pprint("RED","LALALALA")
  for ff in [ff.name for ff in tf.getmembers()]:
    if osp.exists("build/"+ff):
      if osp.isdir("build/"+ff):
        shutil.rmtree("build/"+ff)
      else:
        os.remove("build/"+ff)
  tf.close()
  Logs.pprint("PINK","untar "+what)
  if ctx.exec_command("cd build/;tar -zxf "+what)!=0:
    raise Errors.WafError("Cannot untar "+what)

def installsmthg_post(ctx,where,what,extra_config=""):
  from waflib import Options, Environment,Utils
  CCMACRO = "\"gcc %s\""%ctx.env.mopt
  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
  CPPMACRO = "CPP=\"gcc -E\" CXXCPP=\"g++ -E\" "
  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make;make install"%(where,ctx.env.mprefix,extra_config,CCMACRO, CPPMACRO)
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%what)
  #Logs.pprint("GREEN","You can now run ./waf configure, adding the following option '--%s_islocal'"%what)

