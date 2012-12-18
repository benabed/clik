#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import re
import os.path as osp
import os
import shutil
import clik.hpy as h5py



def main(argv):
  
  lkl = h5py.File(sys.argv[1],"r")
  cls = lkl["clik/check_param"][:]
  if len(argv)==2:
  	clfile = argv+".cls"
  else:
  	clfile = argv[2]
  cls.tofile(clfile,sep=" ")
    
import sys
if __name__=="__main__":
  main(sys.argv)