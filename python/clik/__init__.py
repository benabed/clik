try:
  from lkl import clik
  _lkl_ok = True
except ImportError,e:
  print "Cannot use clik wrapper (cause = '%s')"%e

import re
import numpy as nm

class transformme:
  def __init__(self,tfunc,pf):
    self.tfunc = tfunc
    self.pf =pf
  def __getattr__(self,val):
    return self.tfunc(self.pf[val])

class miniparse(object):
  def __init__(self, pfn):
    print "read parameter file %s"%pfn
    pff =open(pfn)
    txt = "\n".join([to.split("#")[0] for to in pff])+"\n"
    pf = dict(re.findall("(?<!#)(\w+)\s*=\s*(.+?)\n",txt))
    self.pf = pf

  def __contains__(self,val):
    return val in self.pf

  @property
  def int(self):
    return transformme(int,self.pf)

  @property
  def int_array(self):
    return transformme(lambda vl:nm.array([int(v) for v in vl.split() if v]),self.pf)

  @property
  def float(self):
    return transformme(float,self.pf)

  @property
  def float_array(self):
    return transformme(lambda vl:nm.array([float(v) for v in vl.split() if v]),self.pf)

  @property
  def str(self):
    return transformme(lambda x:x,self.pf)

  @property
  def str_array(self):
    return transformme(lambda vl:vl.split(),self.pf)

  def __getattr__(self,val):
    return getattr(self.str,val)