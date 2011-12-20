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
import h5py
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
  lmin_143x143,lmax_143x143,lmin_217x217,lmax_217x217,lmin_143x217,lmax_143x217,np_143x143,np_217x217,np_143x217,nX = fcs.read(" i32"*10)
  X = fcs.read("%df64"%nX)
  fcs.read()
  c_inv = fcs.read("%df64"%(nX)**2)
  
  fsz = forfile(pars.input_szfile)
  lmax_sz = fsz.read('i32')
  sz_temp = fsz.read("%df64"%(lmax_sz+1))

  lmin = nm.min([lmin_143x143,lmin_143x217,lmin_217x217])
  lmax = nm.max([lmax_143x143,lmax_143x217,lmax_217x217])

  hascl = [0]*6
  hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = php.add_lkl_generic(root_grp,"CAMspec",1,hascl,lmax,lmin)

  lkl_grp.attrs["lmin_143x143"] = lmin_143x143
  lkl_grp.attrs["lmin_143x217"] = lmin_143x217
  lkl_grp.attrs["lmin_217x217"] = lmin_217x217
  lkl_grp.attrs["lmax_143x143"] = lmax_143x143
  lkl_grp.attrs["lmax_143x217"] = lmax_143x217
  lkl_grp.attrs["lmax_217x217"] = lmax_217x217
  lkl_grp.attrs["np_143x143"]   = np_143x143
  lkl_grp.attrs["np_143x217"]   = np_143x217
  lkl_grp.attrs["np_217x217"]   = np_217x217

  lkl_grp.attrs["nX"]          = nX
  
  lkl_grp.attrs["lmax_sz"]     = lmax_sz

  lkl_grp.create_dataset('X', data=X)
  lkl_grp.create_dataset('sz_temp', data=sz_temp)
  lkl_grp.create_dataset('c_inv', data=c_inv)
  hf.close()
  
import sys
if __name__=="__main__":
  main(sys.argv)