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
  
  if not Options.options.gcc:
    try:
      ctx.check_tool('icc')
      ctx.check_cc(
        errmsg="failed",msg="Compile a test code with icc",
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
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
  
  if not Options.options.icc:
    try:
      ctx.check_tool('gcc')
      ctx.check_cc(
        errmsg="failed",msg="Compile a test code with gcc",
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
      return
    except:
      if Options.options.gcc:
        raise
      Logs.pprint("PINK", "gcc not found, defaulting to icc")
  ctx.check_tool('icc')
  ctx.check_cc(
        errmsg="failed",msg='Compile a test code with icc',
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')


configure = configure_gccfirst