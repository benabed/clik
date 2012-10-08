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
import clik.smicahlp as smh


  
def main(argv):
  pars = clik.miniparse(argv[1])
  inhf = h5py.File(pars.input_object)
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  inhf.close()

  shutil.copyfile(pars.input_object,pars.res_object)
  outhf = h5py.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]

  crit = pars.criterion
  extra={}
  if crit=="quad":
    if "rq" in pars:
      extra["fid"]=nm.fromfile(pars.rq)
    if "mask" in pars:

      extra["mask"] = nm.asarray(nm.reshape(nm.loadtxt(pars.mask),(lkl_grp.attrs["m_channel_T"],lkl_grp.attrs["m_channel_T"])),nm.int)
      print extra["mask"]
  smh.set_criterion(lkl_grp,crit,**extra)
  outhf.close()  
    
import sys
if __name__=="__main__":
  main(sys.argv)