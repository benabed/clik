#! /usr/bin/env python

import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik

def main(argv):
  lklfile = argv[1]
  lkl = clik.clik(lklfile)
  for clfile in argv[2:]:
    cls = nm.loadtxt(clfile)
    nres = lkl(cls.flat[:])
    print nres

if __name__=="__main__":
  main(sys.argv)