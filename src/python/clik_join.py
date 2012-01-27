#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import clik
import clik.hpy as h5py


def main(argv):
  if len(sys.argv)<4:
    print "usage : %s lkl_file_1 lkl_file_2 [lkl_file_3 ...] result_lkl_file"
    sys.exit(1)
  lkls = [h5py.File(ll)["clik"] for ll in sys.argv[1:-1]]
  reslkl = h5py.File(sys.argv[-1],"w")
  
  nlkl = 0
  lmax = -nm.ones(6)
  resclik = reslkl.create_group("clik")
  name = []
  loc = []
  var = nm.zeros((0,0))

  for lklin in lkls:
    if "clik/prior" in lklin:
      pname = [n.strip() for n in prid.attrs["name"].split()]
      for n in name:
        if n in pname:
          raise Exception("already got a prior on %s"%n)
      ploc = prid["loc"][:]
      pvar = prid["var"][:]
      if len(pvar)==len(ploc):
        pvar = nm.diag(pvar)
      pvar.shape = (len(ploc),-1)
      nvar = nm.zeros((len(loc),len(loc)))
      nvar[:len(loc),:len(loc)] = var
      nvar[len(loc):,len(loc):] = pvar
      var = nvar
      name = list(name) + list(pname)
      loc = nm.concatenate((loc,ploc))


    lmaxin = lklin.attrs["lmax"]
    lmax = nm.max((lmax,lmaxin),0)
    nlklin = lklin.attrs["n_lkl_object"]
    
    for i in range(nlklin):
      grpin = "lkl_%d"%i
      grpout = "lkl_%d"%nlkl
      nlkl+=1
      lklin.copy(lklin[grpin],resclik,grpout)
  
  if len(name):
    prid = resclik.create_group("prior")
    prid.attrs["name"] = php.pack256(*name)
    prid.create_dataset("loc", data=loc.flat[:])
    prid.create_dataset("var", data=var.flat[:])
  
  resclik.attrs["lmax"] = lmax
  resclik.attrs["n_lkl_object"] = nlkl
  
  reslkl.close()
    
if __name__=="__main__":
  main(sys.argv)