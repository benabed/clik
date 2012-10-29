#! $PYTHONEXE
import sys
sys.path = ["$REPLACEPATH"]+sys.path

import h5py 
import clik.smicahlp as hlp
import clik.parobject as php
import numpy as nm
import shutil

def cutlmax(nlmax,infile, outfile):
	shutil.copy(infile,outfile)
	ff = h5py.File(outfile,"r+")
	lmaxs = ff["clik"].attrs["lmax"]
	lmin = ff["clik/lkl_0"].attrs["lmin"]
	blmaxs = ff["clik/lkl_0/bin_lmax"][:]
	idm = blmaxs.searchsorted(nlmax-lmin)
	print "cut at lmax = %d"%(blmaxs[idm]+lmin)
	nnlmax = (blmaxs[idm]+lmin)
	lmaxs[0] = nnlmax
	ff["clik"].attrs["lmax"] = lmaxs
	ff["clik/lkl_0"].attrs["lmax"]=nnlmax
	del ff["clik/lkl_0/bin_lmax"]
	ff["clik/lkl_0/bin_lmax"] = blmaxs[:idm+1]
	blmins = ff["clik/lkl_0/bin_lmin"][:]
	del ff["clik/lkl_0/bin_lmin"]
	ff["clik/lkl_0/bin_lmin"] = blmins[:idm+1]
	wl = ff["clik/lkl_0/bin_ws"][:]
	del ff["clik/lkl_0/bin_ws"]
	ff["clik/lkl_0/bin_ws"] = wl[:nnlmax-lmin+1]
	ff["clik/lkl_0"].attrs["nbins"] = idm+1
	rqh = ff["clik/lkl_0/Rq_hat"][:]
	nch = ff["clik/lkl_0"].attrs["m_channel_T"]
	rqh.shape=[-1,nch,nch]
	rqh = rqh[:idm+1]
	del ff["clik/lkl_0/Rq_hat"]
	ff["clik/lkl_0/Rq_hat"] = rqh.flat[:]
	try :
		blmins = ff["clik/lkl_0/criterion_eig_norm"]
		del ff["clik/lkl_0/criterion_eig_norm"]
		ff["clik/lkl_0/criterion_eig_norm"] = blmins[:idm+1]
	except Exception,e:
		pass
	blmins = ff["clik/lkl_0/wq"][:]
	del ff["clik/lkl_0/wq"]
	ff["clik/lkl_0/wq"] = blmins[:idm+1]

	for i in range(1,ff["clik/lkl_0"].attrs["n_component"]):
		cmpt = "clik/lkl_0/component_%d"%i
		try:
			ff[cmpt].attrs["lmax"]=nnlmax
		except Exception,e:
			pass
		try:
			rqh = ff[cmpt+"/Rq_0"][:]
			rqh.shape=[-1,nch,nch]
			rqh = rqh[:idm+1]
			del ff[cmpt+"/Rq_0"]
			ff[cmpt+"/Rq_0"] = rqh.flat[:]
		except Exception,e:
			pass
		try:
			ng = nm.sum(ff[cmpt].attrs["ngcal"])
			blm = nnlmax-lmin
			if ff[cmpt].attrs["binned"]:
				blm = idm
			tpl = ff[cmpt+"/gcaltpl"][:]
			del ff[cmpt+"/gcaltpl"]
			ff[cmpt+"/gcaltpl"] = tpl[:(blm+1)*ng]
		except Exception,e:
			pass

	try:
		prm = ff["clik/check_param"][:]
		del ff["clik/check_param"]
		del ff["clik"].attrs["check_value"]
		#ff.close()
		#php.add_selfcheck(pars.res_object,cls)
	except Exception,e:
		pass
	ff.close()


import sys
if __name__=="__main__":
  cutlmax(int(sys.argv[1]),sys.argv[2],sys.argv[3])