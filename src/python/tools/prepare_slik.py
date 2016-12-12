#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  
  hascl = [1,1,0,1,0,0]
  hascl = nm.array(hascl,dtype=nm.int)
  lmax = pars.int(default=10).lmax
  lmin = pars.int(default=2).lmin

  
  lkl_grp = php.add_lkl_generic(root_grp,"slik",1,hascl,lmax,lmin)
  
  lkl_grp.attrs["data_lmin"] = 2
  lkl_grp.attrs["data_lmax"] = 12

  lkl_grp.attrs["do_linear"] = pars.int(default=1).do_linear
  lkl_grp.attrs["use_offset_k2"] = pars.int(default=0).do_linear
  lkl_grp.attrs["regcl"] = pars.int(default=1).do_linear

  php.add_pid(lkl_grp,"")

  import tempfile
  dr = tempfile.mkdtemp()
  import os
  print dr
  import shutil
  shutil.copytree(pars.slik_data,dr+"/slik_data")    
  
  php.add_external_data(dr,lkl_grp,tar=True)

  hf.close()

  shutil.rmtree(dr)


import sys
if __name__=="__main__":
  main(sys.argv)