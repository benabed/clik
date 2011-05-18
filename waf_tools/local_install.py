def options(ctx):
  ctx.add_option("--local",action="store_true",default=False,help="install in current directory")

def configure(ctx):
  #install where ?
  ctx.env.mprefix=ctx.env.PREFIX
  import os
  ctx.env.localpref = os.getcwd()
  if ctx.options.local:
    ctx.env.PREFIX=ctx.env.localpref
    ctx.options.prefix = ctx.env.localpref
    ctx.env.mprefix=ctx.env.localpref
    ctx.env.LIBDIR=ctx.env.PREFIX+"/lib"
    ctx.env.BINDIR=ctx.env.PREFIX+"/bin"
