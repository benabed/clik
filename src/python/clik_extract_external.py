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

def main(argv):
  pars = clik.miniparse(argv[1])
  
  inhf = h5py.File(pars.input_object)
  
  dts = inhf["clik/lkl_0/external_data"][:]
  inhf.close()
  if not osp.exists(pars.install_path):
    os.mkdir(pars.install_path)
  f=open(osp.join(pars.install_path,"data.tar"),"w")
  f.write(dts.tostring())
  f.close()
  assert os.system("cd %s;tar xvf data.tar"%pars.install_path)==0
  assert os.system("cd %s;rm -f data.tar"%pars.install_path)==0
  
  shutil.copyfile(pars.input_object,pars.res_object)

  outhf = h5py.File(pars.res_object,"r+")
  del outhf["clik/lkl_0/external_data"]
  outhf["clik/lkl_0"].attrs["external_dir"] = pars.install_path
  outhf.close()

import sys
if __name__=="__main__":
  main(sys.argv)