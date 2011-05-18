def options(ctx):
  ctx.add_option("--gcc",action="store_true",default=False,help="Do not test for icc and only use gcc")

def configure(ctx):
  import Options
  from waflib import Logs
  
  if not Options.options.gcc:
    try:
      ctx.check_tool('icc')
    except:
      Logs.pprint("PINK", "icc not found, defaulting to gcc")
      ctx.check_tool('gcc')
  else:
    ctx.check_tool('gcc')