import os
import os.path as osp
import shutil as shu
import pyfits as pf
import re
import numpy as nm

_metadata = "_mdb"
class File(object):
  def __init__(self,name,mode):
    self._mode = '+'
    if mode=="w":
      self._create(name)
      return
    if mode=="r" or mode=="r-":
      self._name = name
      self._mode="-"
      return
    if mode=="r+":
      self._name=name
    
  def _parsemetadata(self,path=""):
    if not path:
      path = self._name
    f=open(osp.join(path,_metadata))
    dct = {}
    for l in f:
      if not l.strip():
        continue
      id0 = l.find(" ")
      key = l[:id0]
      id1 = l[id0+1:].find(" ") + id0+1
      typ = l[id0+1:id1]
      data = l[id1+1:-1]
      if typ== "int":
        dct[key]=int(data)
        continue
      if typ == "float":
        dct[key] = float(data)
        continue
      if typ == "str":
        dct[key] = data
        continue
      f.close()
      raise TypeError("unknown type '%s' for metadata '%s'"%(typ,key))
    f.close()
    return dct
  def _writemetadata(self,dct,path=""):
    if not path:
      path = self._name
    f=open(osp.join(path,_metadata),"w")
    for k,v in dct.items():
      if type(v)==str:
        typ="str"
        modi = "%s"
      elif type(v) in (bool,int,long):
        typ = "int"
        v = int(v)
        modi = "%d"
      elif type(v) == (float):
        typ="float"
        modi = "%.10g"
      else:
        raise TypeError("bad type")
      f.write(("%s %s "+modi+"\n")%(k,typ,v))
    f.close()

  def remove(self,name):
    if osp.exists(name):
      if osp.isdir(name):
        shu.rmtree(name)
      else:
        os.remove(name)
    else:
      dct = self._parsemetadata(osp.split(name)[0])
      if osp.split(name)[1] in dct.keys():
        del dct[osp.split(name)[1]]
        self._writemetadata(dct,osp.split(name)[0]) 
    
  def _create(self,name):
    if osp.isdir(name):
      shu.rmtree(name)
    os.mkdir(name)
    f=open(osp.join(name,_metadata),"w")
    f.write("")
    f.close()
    self._name = name

  def __getitem__(self,key):
    fkey = osp.join(self._name,key)
    if fkey[-1]=='/':
      fkey = fkey[:-1]
    if osp.exists(fkey):
      if osp.isdir(fkey):
        return File(fkey,"r"+self._mode)
      try:
        return pf.open(fkey)[0].data
      except Exception:
        return open(fkey).read()

    dct = self._parsemetadata(osp.split(fkey)[0])
    return dct[osp.split(fkey)[1]]

  def __setitem__(self,key,value):
    assert self._mode=='+'
    fkey = osp.join(self._name,key)
    if fkey[-1]=='/':
      fkey = fkey[:-1]
    self.remove(fkey)
    if isinstance(value,File):
      
      shu.copytree(value._name,fkey)
      return
    if type(value) in (list,tuple,nm.ndarray):    
      value = nm.array(value)
      if value.dtype==nm.int:
        value = value.astype(nm.long)

      pf.PrimaryHDU(value).writeto(fkey)
      return
    if type(value) == str and ("\n" in value or "\0" in value or len(value)>50):
      f=open(fkey,"w")
      f.write(value)
      f.close()
      return
    dct = self._parsemetadata(osp.split(fkey)[0])
    dct[osp.split(fkey)[1]] = value
    self._writemetadata(dct,osp.split(fkey)[0])    
  
  def create_group(self,name):
    assert self._mode=='+'
    return File(osp.join(self._name,name),"w")
  def create_dataset(self,name,data=None):
    assert data!=None
    self[name] = data
  
  def __delitem__(self,key):
    assert self._mode=='+'
    fkey = osp.join(self._name,key)
    if fkey[-1]=='/':
      fkey = fkey[:-1]

    if osp.exists(fkey):
      self.remove(fkey)
      return 
    dct = self._parsemetadata(osp.split(fkey)[0])
    del dct[osp.split(fkey)[1]]
    self._writemetadata(dct,osp.split(fkey)[0])    

  def copy(self,a,b,c=""):
    if not c:
      self[a] = self[b]
    else:
      b[c]=a
    
  @property
  def attrs(self):        
    return self
  
  def keys(self):
    dct = self._parsemetadata(self._name)
    ls = [el for el in os.listdir(self._name) if el[0]!='.' and el!=_metadata]
    return ls+dct.keys()

  def items(self):
    ks = self.keys()
    return [(k,self[k]) for k in ks]
  

try:
  import h5py
  def hdf2cldf_grp(hdf,fdf):
    # first the metadata
    for kk in hdf.attrs.keys():
      fdf[kk] = hdf.attrs[kk]
    # then the group/data
    for kk in hdf.keys():
      god = hdf[kk]
      if isinstance(god,h5py.Group):
        if not hasattr(fdf,kk):
          fdf.create_group(kk)
        hdf2cldf_grp(god,fdf[kk])
      else:
        fdf[kk] = god[:]

  def hdf2cldf(ffin, ffout):
    hdf = h5py.File(ffin,"r")
    fdf = File(ffout,"w")
    hdf2cldf_grp(hdf,fdf)
except ImportError,e:
  pass
      