_has_h5py = False
try:
	import h5py
	_has_h5py = True
except Exception, e:
	pass
	
import cldf

def File(path,mode="r",ty=None):
	if _has_h5py and type(ty) in dir(h5py):
		try:
			return h5py.File(path,mode)
		except:
			pass
	return cldf.File(path,mode)
