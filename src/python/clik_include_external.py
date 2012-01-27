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
  install_path = inhf["clik/lkl_0"].attrs["external_dir"]
  assert os.system("cd %s;tar cvf data.tar *"%install_path)==0
  f=open(osp.join(install_path,"data.tar"),"r")
  dts = f.read()
  f.close()
  
  shutil.copyfile(pars.input_object,pars.res_object)
  
  outhf = h5py.File(pars.res_object,"r+")
  del outhf["clik/lkl_0"].attrs["external_dir"]
  outhf["clik/lkl_0"].create_dataset("external_data",data=nm.fromstring(dts,dtype=nm.uint8))
  outhf.close()
  

import sys
if __name__=="__main__":
  main(sys.argv)