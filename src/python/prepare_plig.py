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

def expand_bins(bins,ncl):
  rbins = nm.zeros((bins.shape[0]*ncl,bins.shape[1]*ncl))
  for i in range(ncl):
    rbins[bins.shape[0]*i:bins.shape[0]*(i+1),bins.shape[1]*i:bins.shape[1]*(i+1)] = bins
  return rbins

def test_cov_mat_format(fname):
  full = False
  try:
    hdulist = pf.open(fname)
    try:
      dump = hdulist[0].header['DMC_PID']
      full = True
    finally:
      hdulist.close()
  except Exception:
    pass
  return full

def ordering_TEB(nt,np,has_cl):
  m = nt*has_cl[0]+np*has_cl[1]+np*has_cl[2]
  rr=[]
  # TT first
  for m1 in range(nt*has_cl[0]):
    for m2 in range(m1,nt*has_cl[0]):
      rr += [(m1,m2)]
  # EE
  for m1 in range(np*has_cl[1]):
    for m2 in range(m1,np*has_cl[1]):
      rr += [(m1+nt*has_cl[0],m2+nt*has_cl[0])]
  # BB
  for m1 in range(np*has_cl[2]):
    for m2 in range(m1,np*has_cl[2]):
      rr += [(m1+nt*has_cl[0]+np*has_cl[1],m2+nt*has_cl[0]+np*has_cl[1])]
  # TE
  for m1 in range(nt*has_cl[0]):
    for m2 in range(0,np*has_cl[1]):
      rr += [(m1,m2+nt*has_cl[0])]
  # TB
  for m1 in range(nt*has_cl[0]):
    for m2 in range(0,np*has_cl[2]):
      rr += [(m1,m2+nt*has_cl[0]+np*has_cl[1])]
  # EB
  for m1 in range(np*has_cl[1]):
    for m2 in range(0,np*has_cl[2]):
      rr += [(m1+nt*has_cl[0],m2+nt*has_cl[0]+np*has_cl[1])]
  #print nm.array(rr).shape, nt,np,has_cl,rr
  return nm.array(rr)
  
def remove_zero_rowcol(matrix_in, mask):
  idx        = list(set(list(nm.nonzero(mask)[0])))
  matrix_out = nm.zeros([len(idx), len(idx)])
  for i in range(len(idx)):
    matrix_out[i,:] = matrix_in[idx[i],idx]
  return matrix_out

def read_full_cov_mat(fname):
  l_info  = {}
  hdulist = pf.open(fname)
  cov_mat = hdulist['ICOV'].data
  bin_TT  = hdulist['BIN_TT'].data
  bin_EE  = hdulist['BIN_EE'].data
  bin_TE  = hdulist['BIN_TE'].data
  header  = hdulist[0].header
  hdulist.close()
  for key in header.keys():
    if ('LMIN_' in key or 'LMAX_' in key):
      l_info[key] = int(header[key])
  cov_mat /= (1.0E6)**4
  return l_info, bin_TT, bin_EE, bin_TE, cov_mat

def select_channels(l_info, nr_freq, frequencies):
  mask_TP = nm.zeros([2*nr_freq, 2*nr_freq], dtype='int')
  for i, freq1 in enumerate(frequencies):
    for j, freq2 in enumerate(frequencies):
      if i > j:
        continue
      mask_TP[i,j] \
              = (l_info['LMAX_TT_' + freq1 + 'X' + freq2] > 0)
      mask_TP[j,i] = mask_TP[i,j]
      mask_TP[i+nr_freq,j+nr_freq] \
              = (l_info['LMAX_EE_' + freq1 + 'X' + freq2] > 0)
      mask_TP[j+nr_freq,i+nr_freq] = mask_TP[i+nr_freq,j+nr_freq]
      mask_TP[i+nr_freq,j] \
              = (l_info['LMAX_TE_' + freq1 + 'X' + freq2] > 0)
      mask_TP[j+nr_freq,i] = mask_TP[i+nr_freq,j]
      mask_TP[i,j+nr_freq] = mask_TP[i+nr_freq,j]
      mask_TP[j,i+nr_freq] = mask_TP[i+nr_freq,j]
  diag_TT = (sum(mask_TP[:nr_freq,:nr_freq], 0) > 0).astype('int')
  diag_EE = (sum(mask_TP[nr_freq:,nr_freq:], 0) > 0).astype('int')
  diag_TE = (sum(mask_TP[:nr_freq,nr_freq:], 0) > 0).astype('int')
  nT      = sum(diag_TT | diag_TE)
  nP      = sum(diag_EE | diag_TE)
  nTE     = sum(diag_TE)
  has_cl  = [1*(nT > 0), 1*(nP > 0), 0, 1*(nTE > 0), 0, 0]

  frq     = []
  channel = []
  for i, freq in enumerate(frequencies):
    if sum(mask_TP[i,:]) > 0:
      frq.append(float(freq))
      channel.append(freq + 'T')
  for i, freq in enumerate(frequencies):
    if sum(mask_TP[i+nr_freq,:]) > 0:
      frq.append(float(freq))
      channel.append(freq + 'P')
  return nT, nP, has_cl, frq, channel, mask_TP

def get_l_range(l_info, mask_TP, nr_freq, frequencies):
  lmin_TP = -1*nm.ones(mask_TP.shape, dtype='int')
  lmax_TP = -1*nm.ones(mask_TP.shape, dtype='int')
  for i, freq1 in enumerate(frequencies):
    for j, freq2 in enumerate(frequencies):
      if i > j:
        continue
      if mask_TP[i,j] != 0:
        lmin_TP[i,j] = l_info['LMIN_TT_' + freq1 + 'X' + freq2]
        lmin_TP[j,i] = lmin_TP[i,j]
        lmax_TP[i,j] = l_info['LMAX_TT_' + freq1 + 'X' + freq2]
        lmax_TP[j,i] = lmax_TP[i,j]
      if mask_TP[i+nr_freq,j+nr_freq] != 0:
        lmin_TP[i+nr_freq,j+nr_freq] = l_info['LMIN_EE_' + freq1 + 'X' + freq2]
        lmin_TP[j+nr_freq,i+nr_freq] = lmin_TP[i+nr_freq,j+nr_freq]
        lmax_TP[i+nr_freq,j+nr_freq] = l_info['LMAX_EE_' + freq1 + 'X' + freq2]
        lmax_TP[j+nr_freq,i+nr_freq] = lmax_TP[i+nr_freq,j+nr_freq]
      if mask_TP[i+nr_freq,j] != 0:
        lmin_TP[i+nr_freq,j] = l_info['LMIN_TE_' + freq1 + 'X' + freq2]
        lmin_TP[j+nr_freq,i] = lmin_TP[i+nr_freq,j]
        lmin_TP[i,j+nr_freq] = lmin_TP[i+nr_freq,j]
        lmin_TP[j,i+nr_freq] = lmin_TP[i+nr_freq,j]
        lmax_TP[i+nr_freq,j] = l_info['LMAX_TE_' + freq1 + 'X' + freq2]
        lmax_TP[j+nr_freq,i] = lmax_TP[i+nr_freq,j]
        lmax_TP[i,j+nr_freq] = lmax_TP[i+nr_freq,j]
        lmax_TP[j,i+nr_freq] = lmax_TP[i+nr_freq,j]
  try:
    submatrix   = lmin_TP[:nr_freq,:nr_freq]
    min_lmin_TT = min(submatrix[submatrix >= 0])
  except ValueError:
    min_lmin_TT = -1
  submatrix     = lmax_TP[:nr_freq,:nr_freq]
  max_lmax_TT   = max(nm.hstack([-1, submatrix[submatrix >= 0].flatten()]))
  try:
    submatrix   = lmin_TP[nr_freq:,nr_freq:]
    min_lmin_EE = min(submatrix[submatrix >= 0])
  except ValueError:
    min_lmin_EE = -1
  submatrix     = lmax_TP[nr_freq:,nr_freq:]
  max_lmax_EE   = max(nm.hstack([-1, submatrix[submatrix >= 0].flatten()]))
  try:
    submatrix   = lmin_TP[:nr_freq,nr_freq:]
    min_lmin_TE = min(submatrix[submatrix >= 0])
  except ValueError:
    min_lmin_TE = -1
  submatrix     = lmax_TP[:nr_freq,nr_freq:]
  max_lmax_TE   = max(nm.hstack([-1, submatrix[submatrix >= 0].flatten()]))
  lmin          = min(lmin_TP[lmax_TP >= 0])
  lmax          = max(lmax_TP[lmax_TP >= 0])

  l_info['min_lmin_TT'] = min_lmin_TT
  l_info['max_lmax_TT'] = max_lmax_TT
  l_info['min_lmin_EE'] = min_lmin_EE
  l_info['max_lmax_EE'] = max_lmax_EE
  l_info['min_lmin_TE'] = min_lmin_TE
  l_info['max_lmax_TE'] = max_lmax_TE
  l_info['lmin']        = lmin
  l_info['lmax']        = lmax
  return lmin, lmax, lmin_TP, lmax_TP, l_info

def get_l_binning(mask_TP, lmin_TP, lmax_TP, l_info, \
                  bin_TT, bin_TE, bin_EE):
  qmins = nm.zeros(mask_TP.shape, dtype='int')
  qmaxs = nm.zeros(mask_TP.shape, dtype='int')
  if (l_info['min_lmin_TT'] == l_info['lmin']) \
     and (l_info['max_lmax_TT'] == l_info['lmax']):
    bins = bin_TT
  elif (l_info['min_lmin_EE'] == l_info['lmin']) \
     and (l_info['max_lmax_EE'] == l_info['lmax']):
    bins = bin_EE
  elif (l_info['min_lmin_TE'] == l_info['lmin']) \
     and (l_info['max_lmax_TE'] == l_info['lmax']):
    bins = bin_TE
  else:
    print "Error: Using combined binning matrices not implemented"
    quit()
  nr_bins = bins.shape[0]
  lcuts   = nm.zeros(nr_bins+1, dtype='int')
  for i in range(nr_bins):
    lcuts[i] = min(nm.where(bins[i,:] > 0)[0]) + l_info['lmin']
  lcuts[-1] = l_info['lmax'] + 1
  for i in range(mask_TP.shape[0]):
    for j in range(mask_TP.shape[0]):
      if mask_TP[i,j] == 0:
        continue
      qmins[i,j] = nm.where(lcuts == lmin_TP[i,j])[0]
      qmaxs[i,j] = nm.where(lcuts == lmax_TP[i,j]+1)[0]
  qmins = remove_zero_rowcol(qmins, mask_TP)
  qmaxs = remove_zero_rowcol(qmaxs, mask_TP)
  return qmins, qmaxs, nr_bins, bins

def get_power_spectra(fname, nT, nP, nr_freq, mask_TP, l_info, nr_bins, bins):
  cl_raw = read_array(fname)
  try:
    cl_raw.shape = [3001, 2*nr_freq, 2*nr_freq]
  except ValueError:
    print "Error: Power spectrum input file format mismatch"
    quit()
  rqhat     = nm.zeros([nr_bins, nT + nP, nT + nP])
  rqhat_tmp = nm.zeros([nr_bins, mask_TP.shape[0], mask_TP.shape[0]])
  for i in range(mask_TP.shape[0]):
    for j in range(mask_TP.shape[0]):
      rqhat_tmp[:,i,j] = nm.dot(bins, cl_raw[l_info['lmin']:l_info['lmax']+1,i,j])
  for i in range(nr_bins):
    rqhat[i,:,:] = remove_zero_rowcol(rqhat_tmp[i,:,:], mask_TP)
  rqhat *= (1.0E6)**2
  return rqhat

def dump_colorcorrections(color_corr, mask_TP, pars):
  filename  = 'colorcorr_dust.dat'
  write_col = False
  color_in  = pars.str_array.parametric_dot_color
  color_out = ''
  for i, color in enumerate(color_in):
    if color == "%DUST%":
      color_out += ' ' + filename
      write_col = True
    else:
      color_out += ' ' + color
  color_out = color_out.strip()

  if write_col:
    pars.pf['parametric.color'] = color_out
    color_corr_values = []
    for i in range(mask_TP.shape[0]):
      if sum(mask_TP[:,i], 0) > 0:
        color_corr_values.append(color_corr[i])
    handle = open(filename, "w")
    nm.savetxt(handle, color_corr_values, fmt=' %15.7E')
    handle.close()
  return pars

def add_calibration(channel, pars):
  if "calib" in pars and pars.calib.strip():
    return pars
  ref = '143'
  if ((any('T' in entry for entry in channel) and (not ref + 'T' in channel)) or
      (any('P' in entry for entry in channel) and (not ref + 'P' in channel))):
    print "Error: Need {0:3s} GHz channel for calibration".format(ref)
    quit()
  calib_channels = ''
  for freq in channel:
    if ref in freq:
      continue
    if 'P' in freq:
      continue
    calib_channels += ' ' + freq
  calib_channels   = calib_channels.strip()
  if calib_channels:
    pars.pf['calib'] = calib_channels
  return pars

def input_from_cov_mat(pars):
  print "Parsing binning information from covariance matrix"
  frequencies = ['100', '143', '217']
  color_corr  = [1.06881, 1.05195, 1.13962, 1.0, 1.0, 1.0]
  nr_freq     = len(frequencies)

  l_info, bin_TT, bin_EE, bin_TE, cov_mat \
     = read_full_cov_mat(pars.str.mat)

  nT, nP, has_cl, frq, channel, mask_TP \
     = select_channels(l_info, nr_freq, frequencies)

  lmin, lmax, lmin_TP, lmax_TP, l_info \
     = get_l_range(l_info, mask_TP, nr_freq, frequencies)

  qmins, qmaxs, nr_bins, bins \
     = get_l_binning(mask_TP, lmin_TP, lmax_TP, l_info, bin_TT, bin_TE, bin_EE)

  rqhat = get_power_spectra(pars.str.rqhat, nT, nP, nr_freq, mask_TP, \
                            l_info, nr_bins, bins)

  pars = dump_colorcorrections(color_corr, mask_TP, pars)

  pars = add_calibration(channel, pars)

  bins = expand_bins(bins, sum(has_cl))
  Acmb = nm.ones(nT + nP)
  return nT, nP, has_cl, frq, channel, lmin, lmax, nr_bins, bins, \
         qmins, qmaxs, Acmb, rqhat, cov_mat, pars

def input_from_config_file(pars):
  print "Parsing binning information from config file"
  nT = pars.int.nT
  nP = pars.int.nP

  frq = pars.float_array.freq


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

    bins = expand_bins(bins,ncl)

  else:
    lmin = pars.int.lmin
    lmax = pars.int.lmax
    nell = lmax+1-lmin
    nq = lmax+1-lmin
    bins = None

  if "qmins" in pars:
    qmins = pars.int_array.qmin
    qmaxs = pars.int_array.qmax
  else:
    qmins = nm.zeros((len(channel),len(channel)))
    qmaxs = nm.ones((len(channel),len(channel)))*nq

  qmins.shape=((len(channel),len(channel)))
  qmaxs.shape=((len(channel),len(channel)))

  cov_mat = read_array(pars.str.mat)

  rqhat = read_array(pars.str.rqhat)

  Acmb = pars.float_array(default=nm.ones(len(frq))).Acmb
  return nT, nP, has_cl, frq, channel, lmin, lmax, nq, bins, qmins, qmaxs, \
         Acmb, rqhat, cov_mat

def main(argv):
  pars = clik.miniparse(argv[1])

  if test_cov_mat_format(pars.str.mat):
    nT, nP, has_cl, frq, channel, lmin, lmax, nq, bins, qmins, qmaxs, \
        Acmb, rqhat, cov_mat, pars = input_from_cov_mat(pars)
  else:
    nT, nP, has_cl, frq, channel, lmin, lmax, nq, bins, qmins, qmaxs, \
        Acmb, rqhat, cov_mat = input_from_config_file(pars)

  ordering = ordering_TEB(nT,nP,has_cl)

  mask = smh.create_gauss_mask(nq,qmins,qmaxs,nT,nP)

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

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = smh.base_smica(root_grp,has_cl,lmin,lmax,nT,nP,wq,rqhat,Acmb,None,bins)
  smh.set_criterion(lkl_grp,"gauss",mat=cov_mat,mask=mask,ordering=ordering)
  lkl_grp.attrs["dnames"] = php.pack256(*channel)

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
      print pname
      smh.add_parametric_component(lkl_grp,str(pname),frq,keys,lmin,lmax,defaults=defaults,color=colors[ip],rename=rename)

  # Some fix contribution (affected by beam and calib) ?
  if "rq_fix" in pars:
    for rqn in pars.str_array.rq_fix:
      smh.add_cst_component(lkl_grp,read_array(rqn))


  # a gcal component ?
  if "calib" in pars and pars.calib.strip():
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
