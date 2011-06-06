#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik
import h5py 

def main(argv):
  if len(sys.argv)!=2:
    print "usage : %s lkl_file"
    sys.exit(1)
  
  lkl = h5py.File(sys.argv[1])["clik"]
  clikl = clik.clik(sys.argv[1])
  
  print "clik lkl file =  %s"%sys.argv[1]
  print "  number of likelihoods = %d"%lkl.attrs["n_lkl_object"]
  print "  lmax ( "+ " ".join([nl+" = %d"%ll for nl,ll in zip(("TT","EE","BB","TE","TB","EB"),lkl.attrs["lmax"]) if ll >-1])+" )"
  extn = clikl.extra_parameter_names
  print "  number of extra parameters = %d %s"%(len(extn),extn)
  ilkl = 0
  for lkli_n in ("lkl_%d"%v for v in range(lkl.attrs["n_lkl_object"])):
    lkli = lkl[lkli_n]
    print "  %s"%lkli_n
    print "    lkl_type = %s"%lkli.attrs["lkl_type"]
    print "    unit = %g"%lkli.attrs["unit"]
    
    if "lmax" in lkli.attrs:
      lmax = lkli.attrs["lmax"]
      lmin = 0
      if "lmin" in lkli.attrs:
        lmin = lkli.attrs["lmin"]
      ellh = False
    else:
      ell = lkli.attrs["ell"]
      lmax = nm.max(ell)      
      lmin = nm.min(ell)
      ellh = not nm.alltrue((ell[1:]-ell[:-1]) == 1)
    
    print "    "+" ".join([nl+" = [%d , %d]"%(lmin,lmax) for nl,hl in zip(("TT","EE","BB","TE","TB","EB"),lkli.attrs["has_cl"]) if hl ])+" (discontinous)"*ellh
    
    if "wl" in lkli.attrs:
      print "    has window function"
    if "nbins" in lkli.attrs:
      print "    nbins = %d"%lkli.attrs["nbins"]

    extn = clikl.get_extra_parameter_names_by_lkl(ilkl)
    print "    number of extra parameters = %d %s"%(len(extn),extn)
    ilkl +=1
if __name__=="__main__":
  main(sys.argv)