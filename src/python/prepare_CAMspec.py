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

class forfile:
  def __init__(self,fi):
    
    if isinstance(fi,file):
      self.fi = fi
    else :
      self.fi=open(fi)
    self.bf=""
  def read(self,fmt=''):
    if self.bf=='':
      sz = nm.fromstring(self.fi.read(4),dtype=nm.int32)[0]
      #print "want %d bytes"%sz
      self.bf = self.fi.read(sz)
      #print self.bf
      sz2 =nm.fromstring(self.fi.read(4),dtype=nm.int32)[0]
      #print sz2 
      assert sz==sz

    if fmt=='':
      self.bf=''
      return
    res = [self.cvrt(ff) for ff in fmt.strip().split()]
    if len(res)==1:
      return res[0]
    return tuple(res)
  
  def cvrt(self,fmt):
    cmd = re.findall("([0-9]*)([i|f])([0-9]+)",fmt)[0]
    dtype = nm.dtype({"f":"float","i":"int"}[cmd[1]]+cmd[2])
    itm = nm.array(1,dtype=dtype).itemsize
    nelem=1
    if cmd[0]: 
      nelem = int(cmd[0])
    res = nm.fromstring(self.bf[:itm*nelem],dtype=dtype)
    self.bf=self.bf[itm*nelem:]
    if nelem==1:
      return res[0]
    return res
  
  def close(self):
    self.bf=''
    self.fi.close()
  
    

def main(argv):
  pars = clik.miniparse(argv[1])
  fcs = forfile(pars.input_likefile)

  Nspec,nX = fcs.read("i32 i32")

  lminX = []
  lmaxX = []
  np = []
  npt = []
  for i in range(Nspec):
    r = fcs.read("i32 "*4)
    lminX += [r[0]]
    lmaxX += [r[1]]
    np += [r[2]]
    npt += [r[3]]

  X = fcs.read("%df64"%nX)
  fcs.read()
  c_inv = fcs.read("%df64"%(nX)**2)
  
  fsz = forfile(pars.input_szfile_100)
  lmax_sz = fsz.read('i32')
  sz_100 = fsz.read("%df64"%(lmax_sz+1))

  fsz = forfile(pars.input_szfile_143)
  lmax_sz_143 = fsz.read('i32')
  sz_143 = fsz.read("%df64"%(lmax_sz+1))

  assert lmax_sz==lmax_sz_143

  lmin = nm.min(lminX)
  lmax = nm.max(lmaxX)

  hascl = [0]*6
  hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = php.add_lkl_generic(root_grp,"CAMspec",1,hascl,lmax,lmin)

  lkl_grp.attrs["Nspec"] = Nspec
  
  lkl_grp.attrs["lminX"] = lminX
  lkl_grp.attrs["lmaxX"] = lmaxX
  
  lkl_grp.attrs["np"] = np
  lkl_grp.attrs["npt"] = npt

  lkl_grp.attrs["nX"]          = nX
  
  lkl_grp.attrs["lmax_sz"]     = lmax_sz

  lkl_grp.create_dataset('X', data=X)
  lkl_grp.create_dataset('sz_100', data=sz_100)
  lkl_grp.create_dataset('sz_143', data=sz_143)
  lkl_grp.create_dataset('c_inv', data=c_inv)
  hf.close()
  
import sys
if __name__=="__main__":
  main(sys.argv)