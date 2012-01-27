#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import os.path as osp
import os
import shutil
import clik.hpy as h5py
import clik.egfs



def main(argv):
  pars = clik.miniparse(argv[1])
  name = pars.str_array.name
  clikl = clik.clik(pars.input_object)
  extn = clikl.extra_parameter_names
  for n in name:
    if n not in extn:
      raise Exception("extra parameter %s does not exist in likelihood file %s"%(n,pars.input_object))
  shutil.copyfile(pars.input_object,pars.res_object)
  outhf = h5py.File(pars.res_object,"r+")
  loc = pars.float_array.loc
  var = pars.float_array.var
  
  assert len(name)==len(loc)
  assert len(name)==len(var) or len(name)**2==len(var)
  if len(var) == len(loc):
      var = nm.diag(var)
  if "clik/prior" in outhf:
    prid = outhf["clik/prior"]
    pname = [n.strip() for n in prid.attrs["name"].split()]
    for n in name:
      if n in pname:
        raise Exception("already got a prior on %s"%n)
    ploc = prid["loc"][:]
    pvar = prid["var"][:]
    if len(pvar)==len(ploc):
      pvar = nm.diag(pvar)
    pvar.shape = (len(ploc),-1)
    nvar = nm.zeros((len(loc),len(loc)))
    nvar[:len(loc),:len(loc)] = var
    nvar[len(loc):,len(loc):] = pvar
    var = nvar
    name = list(name) + list(pname)
    print name
    loc = nm.concatenate((loc,ploc))
  else:
    prid = outhf.create_group("clik/prior")
  
  var.shape = (len(loc),-1)
  if nm.alltrue(var==nm.diag(nm.diagonal(var))):
    var = nm.diagonal(var)
  prid.attrs["name"] = php.pack256(*name)
  prid.create_dataset("loc", data=loc.flat[:])
  prid.create_dataset("var", data=var.flat[:])
    
  
  if "check_param" in outhf["clik"]:
    cls = outhf["clik/check_param"]
    del(outhf["clik/check_param"])
    del(outhf["clik/check_value"])
    outhf.close()
    php.add_selfcheck(pars.res_object,cls)
 

    
import sys
if __name__=="__main__":
  main(sys.argv)