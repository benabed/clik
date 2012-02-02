def options(ctx):
  import optparse
  grp = ctx.parser.get_option_group("--gfortran")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"compiler options")
  grp.add_option("--gcc",action="store_true",default=False,help="Do not test for icc and only use gcc")
  grp.add_option("--icc",action="store_true",default=False,help="Do not test for gcc and only use icc")
  ctx.add_option_group(grp)  

def configure_iccfirst(ctx):
  import Options
  from waflib import Logs
  ctx.env.has_icc = False
  if not Options.options.gcc:
    try:
      ctx.check_tool('icc')
      ctx.check_cc(
        errmsg="failed",msg="Compile a test code with icc",
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
      ctx.env.has_icc = True
      return
    except:
      if Options.options.icc:
        raise
      Logs.pprint("PINK", "icc not found, defaulting to gcc")
  ctx.check_tool('gcc')
  ctx.check_cc(
        errmsg="failed",msg='Compile a test code with gcc',
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')


def configure_gccfirst(ctx):
  import Options
  from waflib import Logs
  from waflib import  Context
  from waflib import Errors

  ctx.env.has_icc = False
  import re  
  if not Options.options.icc:
    try:
      ctx.check_tool('gcc')
      ctx.start_msg("Check gcc version") 
      v90 = ctx.cmd_and_log(ctx.env.CC[0]+" --version",quiet=Context.STDOUT).split("\n")[0].strip()
      version90 = re.findall("(4\.[0-9]\.[0-9])",v90)
      if len(version90)<1:
        #Logs.pprint("PINK","Can't get gfortran version... Let's hope for the best")
        ctx.end_msg("not found, let's hope for the best...",color="PINK")
      else:
        version90 = version90[0]
        vmid = int(version90.split(".")[1])
        if vmid<2:
          ctx.end_msg(v90,color="YELLOW")
          raise Errors.WafError("gcc version need to be above 4.2 got %s"%version90)
        ctx.end_msg(v90)
      ctx.check_cc(
        errmsg="failed",msg="Compile a test code with gcc",
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
      return
    except Exception,e:
      if Options.options.gcc:
        raise
      Logs.pprint("PINK", "gcc not found, defaulting to icc (cause : %s)"%e)
  ctx.check_tool('icc')
  ctx.check_cc(
        errmsg="failed",msg='Compile a test code with icc',
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
  ctx.env.has_icc = True

configure = configure_gccfirst