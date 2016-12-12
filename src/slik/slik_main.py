#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path
import numpy as nm
import sys
import os

import slik
import clik.hpy as hpy

def main(argv):
  
  #print >>sys.stderr,"BEFORE"
  #print >>sys.stderr,"ivi"
  lklfile = sys.stdin.readline().strip()
  #print >>sys.stderr,"lklfile"
  lkl = hpy.File(lklfile)["lkl_0"]

  slk = slik.stevenlike(lmin=lkl["data_lmin"],lmax=lkl["data_lmax"],loadroot = lkl._name+"/_external/slik_data/", 
                     do_linear=lkl["do_linear"]==1,
                     use_offset_k2=lkl["use_offset_k2"]==1,
                     lminlike=lkl["lmin"],
                     lmaxlike=lkl["lmax"],
                     regcl=lkl["regcl"]==1)

  bl = slk.b_l
  lmin = slk.lmin
  lmax = slk.lmax
  lminlike = slk.lminlike
  lmaxlike = slk.lmaxlike

  ltot = 3*(lmaxlike-lminlike+1)  

  pls = slk.fidcls.copy()
  #print >>sys.stderr,"READY"
  
  print "rep: READY"
  sys.stdout.flush()

  cnt = 0
  rls=[]
  while(True):
    vo = sys.stdin.readline()
    if vo.strip()=="stop":
      #print>>sys.stderr, "bye"
      break
    if vo.strip()[:4]!="val:":
      continue
    else:
      rls+= [float(vo.strip()[5:])]
      #print "got %g"%rls[-1]
      sys.stdout.flush()
      cnt +=1
    if cnt!=ltot:
      continue
    #print "got them all (%d values)"%cnt
    sys.stdout.flush()
    cls = rls
    rls = []
    cnt = 0
    pos = 0
    for i in [0,2,1]:
      for j in range(lmaxlike+1-lminlike):
        pls[(j+lmin-lminlike)*3+i] = cls[pos]*bl[j+lminlike]**2
        pos+=1
    
    #print "reorg now compute"
    sys.stdout.flush()
    lkl = slk.lkl(pls)
    #print "lkl is %g"%lkl
    print "rep: READY"
    sys.stdout.flush()
    print lkl
    sys.stdout.flush()
    




import sys
if __name__=="__main__":
 main(sys.argv)