import os.path as osp

if osp.exists(osp.join(osp.dirname(__file__),"lkl.pyx")):
  raise ImportError("Cannot import clik python wrapper from the source directory.\nMake sure that you have compiled and installed clik and then\nrun python from another directory.")

try:
  from lkl import clik
  _lkl_ok = True
except ImportError,e:
  print "Cannot use clik wrapper (cause = '%s')"%e

import re
import numpy as nm

from miniparse import miniparse

