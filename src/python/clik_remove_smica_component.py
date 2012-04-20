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
import clik.smicahlp as slp

def main(argv):
  pars = clik.miniparse(argv[1])
  
  inhf = h5py.File(pars.input_object)
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  inhf.close()
  
  shutil.copyfile(pars.input_object,pars.res_object)
  outhf = h5py.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]
  del lkl_grp["component_%d"%component_index]
  
  php.remove_selfcheck(root_grp = outhf["clik"])
  outhf.close()

import sys
if __name__=="__main__":
  main(sys.argv)