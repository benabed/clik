try:
  from lkl import clik
  _lkl_ok = True
except ImportError,e:
  print "Cannot use clik wrapper (cause = '%s')"%e
