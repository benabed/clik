def options(ctx):
  ctx.add_option("--local",action="store_true",default=False,help="install in current directory")

def configure(ctx):
  #install where ?
  ctx.env.mprefix=ctx.env.PREFIX
  import os
  import os.path as osp
  ctx.env.localpref = os.getcwd()
  if ctx.options.local:
    ctx.env.PREFIX=ctx.env.localpref
    ctx.options.prefix = ctx.env.localpref
    ctx.env.mprefix=ctx.env.localpref
    ctx.env.LIBDIR=osp.join(ctx.env.PREFIX,"lib")
    ctx.env.BINDIR=osp.join(ctx.env.PREFIX,"bin")
  
  if not os.path.exists(ctx.env.LIBDIR):
    os.mkdir(ctx.env.LIBDIR)
    
  if not os.path.exists(ctx.env.BINDIR):
    os.mkdir(ctx.env.BINDIR)
    
  if not os.path.exists(osp.join(ctx.env.PREFIX,"include")):
    os.mkdir(osp.join(ctx.env.PREFIX,"include"))