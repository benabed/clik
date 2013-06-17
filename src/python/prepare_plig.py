#! PYTHONEXE

import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import clik.hpy as h5py
import clik.smicahlp as smh

    
def read_array(fname):
  try:
    return pf.open(fname)[0].data
  except Exception:
    return nm.loadtxt(fname)

def main(argv):
  pars = clik.miniparse(argv[1])

  frq = pars.float_array.freq
  channel = frq

  lmin = pars.int.lmin
  lmax = pars.int.lmax

  nq = lmax+1-lmin
  nell = nq

  nT = pars.int.nT
  nP = pars.int.nP

  has_cl = pars.int_array.has_cl

  ncls = nm.sum(has_cl)

  
  if "bins.limit" in pars:
    blims = pars.int_array.bins_dot_limit
    nq = len(blims)-1
    qwgh = pars.float_array.bins_dot_weight
    bins = nm.zeros(nq,nell)
    bm = blims[0]
    wi = 0
    lmin = pars.int(default=blims[0]).lmin
    lmax = pars.int(default=blims[-1]-1).lmax

    for i,bM in enumerate(blims[1:]):
      nb = bM - bm
      bins[i,bm-lmin:bM-lmin] = qwgh[wi:wi+nb]
      wi+=nb
    bins = nm.tile(bins,(ncl,ncl))
  else:
    lmin = pars.int.lmin
    lmax = pars.int.lmax
    nq = lmax+1-lmin
    bins = None

  qmins = pars.int_array.qmin
  qmaxs = pars.int_array.qmax

  qmins.shape=((len(channel),len(channel)))
  qmaxs.shape=((len(channel),len(channel)))

  mask = smh.create_gauss_mask(nq,qmins,qmaxs)

  wq = nm.ones(nq) *1.


  #if "rqhat" in pars:
  #  rqhat = read_array(pars.str.rqhat)
  #else:
  #  ordering_cl = [[int(v) for v in l.split("x")] for l in pars.str_array.cl_file_order]
  #  rqhat = nm.zeros((nq,(len(channel),len(channel))))
  #  cls = [read_cl(clc,qmins[o[0],o[1]],qmaxs[o[0],o[1]]) for o, clc in zip(ordering_cl,pars.str_array.cl_file)]
  #  for cl,o in zip(pars.str_array.cl_file,ordering_cl):
  #    cls = read_cl(cl,qmins[o[0],o[1]],qmaxs[o[0],o[1]])
  #    rqhat[qmins[o[0],o[1]]:qmaxs[o[0],o[1]],o[0],o[1]] = cls
  #    rqhat[qmins[o[0],o[1]]:qmaxs[o[0],o[1]],o[1],o[0]] = cls

  rqhat = read_array(pars.str.rqhat)

  Acmb = pars.float_array(default=nm.ones(len(channel))).Acmb

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = smh.base_smica(root_grp,has_cl,lmin,lmax,nT,nP,wq,rqhat,Acmb,None,bins)
  smh.set_criterion(lkl_grp,"gauss",mat=read_array(pars.str.mat),mask=mask)


  # Some noise ?
  if "rq_noise" in pars:
    for rqn in pars.str_array.rq_noise:
      smh.add_cst_component(lkl_grp,read_array(rqn))

  # a gcal component ?
  #TBD

  # parametric components ?
  if "parametric" in pars:
    defaults = {}
    if "parametric.default.parameter" in pars:
      defaults = dict(zip(pars.str_array.parameteric_dot_default_dot_parameter,pars.str_array.parameteric_dot_default_dot_value))
    rename = {}
    if "parametric.rename_dot_from" in pars:
      rename = dict(zip(pars.str_array.parameteric_dot_rename_dot_from,pars.str_array.parameteric_dot_rename_dot_to))
    
    keys = pars.parameteric_dot_parameter
    colors = [none]*1000
    if "parametric.color" in pars:
      colors = []
      for cl in pars.str_array.parameteric_dot_color:
        if cl.lower()=="none":
          colors += [nm.ones(len(channel))]
        else:
          colors += [read_array(cl)]
    for ip,pname in enumerate(pars.str_array.parametric):
      smh.add_parametric_component(lkl_grp,pname,frq,keys,lmin,lmax,defaults=defaults,color=colors[i],rename=rename)


  hf.close()
  





import sys
if __name__=="__main__":
  main(sys.argv)