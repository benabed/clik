#! /usr/bin/env python

import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik
import h5py 

def main(argv):
  if len(sys.argv)<4:
    print "usage : %s lkl_file_1 lkl_file_2 [lkl_file_3 ...] result_lkl_file"
    sys.exit(1)
  lkls = [h5py.File(ll)["clik"] for ll in sys.argv[1:-1]]
  reslkl = h5py.File(sys.argv[-1],"w")
  
  nlkl = 0
  lmax = -nm.ones(6)
  resclik = reslkl.create_group("clik")

  for lklin in lkls:
    lmaxin = lklin.attrs["lmax"]
    lmax = nm.max((lmax,lmaxin),0)
    nlklin = lklin.attrs["n_lkl_object"]
    
    for i in range(nlklin):
      grpin = "lkl_%d"%i
      grpout = "lkl_%d"%nlkl
      nlkl+=1
      lklin.copy(lklin[grpin],resclik,grpout)
  
  resclik.attrs["lmax"] = lmax
  resclik.attrs["n_lkl_object"] = nlkl
  
  reslkl.close()
    
if __name__=="__main__":
  main(sys.argv)