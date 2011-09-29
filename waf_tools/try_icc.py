def options(ctx):
  import optparse
  grp = ctx.parser.get_option_group("--gfortran")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"compiler options")
  grp.add_option("--gcc",action="store_true",default=False,help="Do not test for icc and only use gcc")
  ctx.add_option_group(grp)  

def configure(ctx):
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
      Logs.pprint("PINK", "icc not found, defaulting to gcc")
  ctx.check_tool('gcc')
  ctx.check_cc(
        errmsg="failed",msg='Compile a test code with gcc',
        mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')