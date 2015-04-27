#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import os.path as osp
import os
import shutil
import clik.hpy as hpy


def change_gibbs_gauss(inhf,lklfile,outfile,lmin,lmax):
  if lmin == -1 or lmin>lmax:
    lmin = inhf["clik/lkl_0/lmin"]
  if lmax == -1:
    lmax = inhf["clik/lkl_0/lmax"]
  if lmax>249:
    print "not possible"
    sys.exit(-1)
  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  outhf["clik/lmax"] = [lmax,-1,-1,-1,-1,-1]
  outhf["clik/lkl_0/lmin"] = lmin
  outhf["clik/lkl_0/lmax"] = lmax
  outhf["clik/lkl_0/delta_l"] = lmax
  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()

def change_smica(inhf,outfile,lmin,lmax):
  pass

def main(argv):

  lklfile = sys.argv[1]
  lmin = int(sys.argv[2])
  lmax = int(sys.argv[3])
  outfile = sys.argv[4]

  inhf = hpy.File(lklfile)
  ty = inhf["clik/lkl_0/lkl_type"]
  if ty not in ("smica","gibbs_gauss"):
    print "can only change lmin and lmax for plik and commander TT likelihoods"
    sys.exit(-1)
  if ty =="smica":
    change_smica(inhf,lklfile,outfile,lmin,lmax)
  else:
    change_gibbs_gauss(inhf,lklfile,outfile,lmin,lmax) 

    
import sys
if __name__=="__main__":
  main(sys.argv)