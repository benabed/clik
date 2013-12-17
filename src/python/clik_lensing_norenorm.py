#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import clik.hpy as hpy
import os.path as osp

def main(argv):

  lkl = clik.clik_lensing(argv[1])
  del(lkl)
  rr = open(argv[1]).read()
  if ("# format: mono") in rr:
    itype = 1
  else:
    itype = 2
  root = hpy.File(argv[2],"w")
  root.create_group("clik_lensing")
  root["clik_lensing/renorm"]=0
  root["clik_lensing/itype"]=itype
  hpy.copyfile(argv[1],argv[2]+"/clik_lensing/"+osp.basename(argv[1]))
  root["clik_lensing/filename"]=osp.basename(argv[1])
    
import sys
if __name__=="__main__":
  main(sys.argv)