#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik
import clik.hpy as h5py

def main(argv):
  if len(sys.argv)!=2:
    print "usage : %s lkl_file"
    sys.exit(1)
  
  clikl = clik.clik(sys.argv[1])
  extn = clikl.extra_parameter_names
  
  lkl = h5py.File(sys.argv[1],"r")["clik"]
  
  print "clik lkl file =  %s"%sys.argv[1]
  print "  number of likelihoods = %d"%lkl.attrs["n_lkl_object"]
  print "  lmax ( "+ " ".join([nl+" = %d"%ll for nl,ll in zip(("TT","EE","BB","TE","TB","EB"),lkl.attrs["lmax"]) if ll >-1])+" )"
  print "  number of varying extra parameters %d"%(len(extn))
  for n in extn:
    print "    %s"%n
  if "prior" in lkl:
    print "  gaussian priors on %s"%lkl["prior"].attrs["name"]
    loc = lkl["prior/loc"][:]
    print "  at \n    %s"%" ".join([str(l) for l in loc])
    var = lkl["prior/var"][:]
    if len(var)==len(loc):
      var = nm.diag(var)
    var.shape=(len(loc),-1)
    print "  with variance"
    print "\n".join(["    "+" ".join([str(v) for v in vl]) for vl in var])
  if "default" in lkl:
    loc = lkl["default/loc"][:] 
    print "  number of fixed parameters = %d"%len(loc)
    nms = lkl["default"].attrs["name"]
    nms = [nms[i*256:i*256+256].strip() for i in range(len(loc))]
    for n,l in zip(nms,loc):
      print "    %s = %g"%(n,l)

  ilkl = 0
  for lkli_n in ("lkl_%d"%v for v in range(lkl.attrs["n_lkl_object"])):
    lkli = lkl[lkli_n]
    print "\n  %s"%lkli_n
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
    if lkli.attrs["lkl_type"]=="smica":
      print "    component 0 : CMB"
      for nc in range(1,lkli.attrs["n_component"]):
        print "    component %d : %s"%(nc,lkli["component_%d"%nc].attrs["component_type"])

    extn = clikl.get_extra_parameter_names_by_lkl(ilkl)
    print "    number of extra parameters = %d %s"%(len(extn),extn)
    ilkl +=1
if __name__=="__main__":
  main(sys.argv)