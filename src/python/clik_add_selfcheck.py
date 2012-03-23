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
import clik.egfs



def main(argv):
  pars = clik.miniparse(argv[1])
  cls = nm.loadtxt(pars.input_cl)
  res = php.add_selfcheck(pars.input_object,cls)
  print "lkl for init cl %g"%res 
    
import sys
if __name__=="__main__":
  main(sys.argv)