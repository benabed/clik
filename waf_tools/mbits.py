def options(ctx):
  ctx.add_option("--m32",action="store_true",default=False,help="compile & link in 32bits")
  ctx.add_option("--m64",action="store_true",default=False,help="compile & link in 64bits")
  
def configure(ctx):
  import sys
  #32 bits
  if sys.platform.lower()=="darwin":
    mopt = ""
    if ctx.options.m64:
      mopt += "-arch x86_64 "
    if ctx.options.m32:    
      mopt += "-arch i386 "    
  else:
    mopt = "-m64"
    if ctx.options.m32:
      mopt = "-m32"
      
  ctx.env.mopt=mopt
  ctx.env.append_value('CCFLAGS',mopt.split())
  ctx.env.append_value('LINKFLAGS',mopt.split())
