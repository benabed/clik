#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re


def exists(partxt,data_dir,ENTRY):
  ent = re.findall(ENTRY+"\s*=\s*(.+)",partxt)[0]
  return osp.exists(osp.join(data_dir,ent))
  
def main(argv):
  pars = clik.miniparse(argv[1])

  bopix_dir = osp.realpath(pars.bopix_data)
  partxt=open(osp.join(bopix_dir,"input_file.txt")).read()

  # get data directory
  data_dir = osp.join(bopix_dir,re.findall("BOPIX_DATA_DIR\s*=\s*(.+)",partxt)[0])
  assert osp.exists(data_dir)
  
  for fi in ["BOPIX_CL_FILENAME","COV_FILE","WINDOW_FILE"]
    assert exists(partxt,data_dir,fi)
  
  if not exists(partxt,data_dir,"MAP_FILE"):
    assert exists(partxt,data_dir,"MASKED_MAP_P")
    assert exists(partxt,data_dir,"MASKED_MAP_T")

  if not exists(partxt,data_dir,"MASKFILE"):
    assert exists(partxt,data_dir,"MASKFILE_P")
    assert exists(partxt,data_dir,"MASKFILE_T")

  assert os.system("tar cvf /tmp/bopix_data.tar %s"%bopix_dir)==0
  dts = open("/tmp/bopix_data.tar").read()
  assert os.system("rm /tmp/bopix_data.tar")==0
  
  lmax = int(re.findall("BOPIX_CL_LMAX\s*=\s*(.+)",partxt)[0])
  lmin = 0

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = nm.array([1,1,0,1,0,0],dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"bopix",1,hascl,lmax,lmin)
  lkl_grp.create_dataset("external_data",data=nm.frombuffer(dts,dtype=nm.uint8))
  
  hf.close()
  
  
import sys
if __name__=="__main__":
  main(sys.argv)