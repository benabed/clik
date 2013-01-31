_has_h5py = False
try:
	import h5py
	_has_h5py = True
except Exception, e:
	#print e
	pass

_has_cldf = False
try:
	import cldf
	_has_cldf = True
except Exception,e:
	pass

import os.path as osp

def is_h5py_object(oo):
	if not _has_h5py:
		return False
	if not _has_cldf:
		return True
	return not(type(oo) == cldf.File)

def File(path,mode="r",ty=None):
	if _has_h5py and (is_h5py_object(ty) or ty==None):
		try:
			if osp.exists(path) and osp.isdir(path):
				pass
			else:
				return h5py.File(path,mode)
		except Exception,e:
			#print e
	return cldf.File(path,mode)
