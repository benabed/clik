import autoinstall_lib as atl

version = "hdf5-1.8.6"
tool = "hdf5"

#print "-> loading %s autoinstall (using version %s)"%(tool,version)

def options(opt):
  atl.add_lib_option(tool,opt,install=True)
  
def configure(ctx):
  atl.conf_lib(ctx,tool,["hdf5","hdf5_hl"],"H5Fcreate","hdf5.h",defines="HAS_HDF5",install=installhdf5)  

def installhdf5(ctx):
  filen = version+".tar.gz"
  atl.installsmthg_pre(ctx,"http://www.hdfgroup.org/ftp/HDF5/prev-releases/"+version+"/src/"+filen,filen)
  atl.installsmthg_post(ctx,version,tool)
