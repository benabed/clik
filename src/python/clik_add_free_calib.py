#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import clik
import clik.hpy as hpy



def main(argv):
  pars = clik.miniparse(argv[1])
  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
  lkli = pars.int(default=0).lkl_id
  outhf["clik/lkl_%d"%lkli].attrs["free_calib"] = pars.str.parname
  if "check_param" in outhf["clik"]:
    del(outhf["clik/check_param"])
    del(outhf["clik/check_value"])
    
  outhf.close()

    
import sys
if __name__=="__main__":
  main(sys.argv)