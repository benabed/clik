import numpy as nm
import h5py 



def baseCreateParobject(parobject):
  # init file
  hf = h5py.File(parobject, 'w')
  root_grp = hf.create_group("clik")
  
  # fill general info
  root_grp.attrs["n_lkl_object"] = 0
  root_grp.attrs["lmax"] = [-1,-1,-1,-1,-1,-1]
  
  return root_grp,hf
  
def add_external_data(directory,lkl_grp,tar=False):
  import os.path as osp
  import tempfile
  import tarfile
  import numpy as nm
  if not tar:
    lkl_grp.attrs["external_dir"] = osp.realpath(directory)
  else:
    tmp = tempfile.TemporaryFile()
    tartmp = tarfile.TarFile(mode = "w", fileobj=tmp)
    tartmp.add(directory)
    tartmp.close()
    tmp.seek(0)
    dat = nm.frombuffer(tmp.read(),dtype=nm.uint8)
    lkl_grp.create_dataset("external_data",data=dat.flat[:])
    tmp.close()
    
    
def add_lkl_generic(root_grp,lkl_type,unit,has_cl,lmax=-1,lmin=-1,ell=None,wl=None,nbins=0,bins=None,compress_bns=True):
  ilkl = root_grp.attrs["n_lkl_object"]
  lmaxs = root_grp.attrs["lmax"]
  name = "lkl_%d"%ilkl
  
  lkl_grp = root_grp.create_group(name)
  lkl_grp.attrs["lkl_type"] = lkl_type
  lkl_grp.attrs["unit"]     = unit
  lkl_grp.attrs["has_cl"]   = has_cl
  ncl = nm.sum(has_cl)
  
  assert not (lmax == -1 and ell == None)
  assert not (lmax != -1 and ell != None)

  if ell != None:
    lkl_grp.attrs["ell"] = nm.sort(ell)
    lmax = max(ell)
  
  if lmax>-1:
    lkl_grp.attrs["lmax"] = lmax
    ell = nm.arange(lmax+1)
    if lmin>-1:
      lkl_grp.attrs["lmin"] = lmin
      ell = ell[lmin:]
      
  
  
  if nbins>0:
    lkl_grp.attrs["nbins"] = int(nbins)
    if compress_bns==True:
      b_ws,blmin,blmax = compress_bins(bins)
      if b_ws.size+2*blmin.size<bins.size:
        print "compressing bins"
        lkl_grp.create_dataset("bin_ws",data=b_ws.flat[:])
        lkl_grp.create_dataset("bin_lmin",data=blmin.flat[:])
        lkl_grp.create_dataset("bin_lmax",data=blmax.flat[:])
      else:
        compress_bns==False
    if compress_bns==False:
      lkl_grp.create_dataset("bins",data=bins.flat[:])

  if wl!=None:
    lkl_grp.attrs["wl"] = wl
    
  ilkl+=1
  root_grp.attrs["n_lkl_object"] = ilkl
  lmaxs = [max(lm,((lmax+1)*hcl)-1) for lm,hcl in zip(lmaxs,has_cl)]
  root_grp.attrs["lmax"] = lmaxs
    
  return lkl_grp
  
def compress_bins(bins):
  mins = bins!=0
  l = nm.arange(bins.shape[-1])
  blmin,blmax = nm.array([(l[ins][0],l[ins][-1]) for ins in mins]).T
  b_ws = bins[mins]
  return b_ws,blmin,blmax
  
def uncompress_bins(shape,b_ws,blmin,blmax):
  bins = nm.zeros(shape)
  lc = 0
  for i in range(shape[0]):
    bsz = blmax[i]-blmin[i]+1
    bins[i,blmin[i]:blmax[i]+1] = b_ws[lc:lc+bsz]
    lc+=bsz
  return bins