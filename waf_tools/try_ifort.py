import sys
def options(ctx):
  ctx.add_option("--gfortran",action="store_true",default=False,help="Do not test for ifort and only use gfortran")
  ctx.add_option("--fortran_flagline",action="store",default="",help="flagline to link fortran object using ld")

def configure(ctx):
  import Options
  from waflib import Logs
  ifort = False
  
  if not Options.options.gcc:
    try:
      ctx.check_tool('ifort')
      if sys.platform.lower()=="darwin":
        ctx.env.LINKFLAGS_fcshlib = ['-dynamiclib']
      ifort = True
    except:
      Logs.pprint("PINK", "icc not found, defaulting to gcc")
      ctx.check_tool('gfortran')
  else:
    ctx.check_tool('gfortran')
    
  if sys.platform.lower()=="darwin":
    ctx.env.fcshlib_PATTERN = 'lib%s.dylib'

  if ctx.options.fortran_flagline:
    conf.parse_flags(ctx.options.fortran_flagline,uselib="fc_runtime")
  else:
    if ifort:
      import os.path as osp
      ifort_path = osp.dirname(osp.realpath(ctx.env.FC))
      #print ifort_path
      try:
        f=open(osp.join(ifort_path,'ifortvars_intel64.sh'))
      except:
        f=open(osp.join(ifort_path,'ifortvars_ia32.sh'))
      txt = f.read()
      f.close()
      import re
      #print txt
      res = re.findall("DYLD_LIBRARY_PATH\s*=\s*\"(.+)\"",txt)[0]
      for pth in res.split(":"):
        ctx.env.append_value("LIBPATH_fc_runtime",pth)
        ctx.env.append_value("RPATH_fc_runtime",pth)
      ctx.env.append_value("LIB_fc_runtime",["ifcore","intlc","ifport","imf","irc","svml","iomp5"])
      