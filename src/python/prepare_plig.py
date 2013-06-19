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
import pyfits as pf
    
def read_array(fname):
  try:
    pfits = pf.open(fname)
    ii=0
    while pfits[ii].data == None:
      ii+=1
    return pfits[ii].data
  except Exception:
    return nm.loadtxt(fname)

def main(argv):
  pars = clik.miniparse(argv[1])

  frq = pars.float_array.freq

  
  nT = pars.int.nT
  nP = pars.int.nP

  channel = pars.str_array.channel
  has_cl = pars.int_array.has_cl

  ncl = nm.sum(has_cl)

  
  if "bins.limit" in pars:
    blims = pars.int_array.bins_dot_limit
    lmin = pars.int(default=blims[0]).lmin
    lmax = pars.int(default=blims[-1]-1).lmax
    nell = lmax+1-lmin
    nq = len(blims)-1
    qwgh = pars.float_array.bins_dot_weights
    bins = nm.zeros((nq,nell),dtype=nm.double)
    bm = blims[0]
    wi = 0
    lmin = pars.int(default=blims[0]).lmin
    lmax = pars.int(default=blims[-1]-1).lmax

    for i,bM in enumerate(blims[1:]):
      nb = bM - bm
      bins[i,bm-lmin:bM-lmin] = qwgh[wi:wi+nb]
      wi+=nb
      bm=bM
    bins = nm.tile(bins,(ncl,1))
  else:
    lmin = pars.int.lmin
    lmax = pars.int.lmax
    nell = lmax+1-lmin
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

  Acmb = pars.float_array(default=nm.ones(len(frq))).Acmb

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = smh.base_smica(root_grp,has_cl,lmin,lmax,nT,nP,wq,rqhat,Acmb,None,bins)
  smh.set_criterion(lkl_grp,"gauss",mat=read_array(pars.str.mat),mask=mask)
  lkl_grp.attrs["dnames"] =  php.pack256(*channel) 

  # parametric components ?
  if "parametric" in pars:
    defaults = {}
    if "parametric.default.parameter" in pars:
      defaults = dict(zip(pars.str_array.parametric_dot_default_dot_parameter,pars.str_array.parametric_dot_default_dot_value))
    rename = {}
    if "parametric.rename.from" in pars:
      rename = dict(zip(pars.str_array.parametric_dot_rename_dot_from,pars.str_array.parametric_dot_rename_dot_to))
    
    keys = pars.str_array.parametric_dot_parameter
    colors = [None]*1000
    if "parametric.color" in pars:
      colors = []
      for cl in pars.str_array.parametric_dot_color:
        if cl.lower()=="none":
          colors += [nm.ones(len(frq))]
        else:
          colors += [read_array(cl)]
    for ip,pname in enumerate(pars.str_array.parametric):
      smh.add_parametric_component(lkl_grp,str(pname),frq,keys,lmin,lmax,defaults=defaults,color=colors[ip],rename=rename)

  # Some fix contribution (affected by beam and calib) ?
  if "rq_fix" in pars:
    for rqn in pars.str_array.rq_fix:
      smh.add_cst_component(lkl_grp,read_array(rqn))


  # a gcal component ?
  if "calib" in pars:
    names = ["calib_"+v for v in pars.str_array.calib]
    smh.add_calTP_component(lkl_grp,names)
  if "beammode.select" in pars:
    names = ["beammode_"+v for v in pars.str_array.beammode_dot_select]
    tpl = [read_array(v) for v in pars.str_array.beammode_dot_data]
    smh.add_gcal2_component(lkl_grp,names,tpl)

  # Some noise ?
  if "rq_noise" in pars:
    for rqn in pars.str_array.rq_noise:
      smh.add_cst_component(lkl_grp,read_array(rqn))
  
  hf.close()
  



import sys
if __name__=="__main__":
  main(sys.argv)