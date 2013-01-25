_has_h5py = False
try:
	import h5py
	_has_h5py = True
except Exception, e:
	print e
	pass
	
import cldf

def File(path,mode="r",ty=None):
	if _has_h5py and (type(ty) in dir(h5py) or ty==None):
		try:
			return h5py.File(path,mode)
		except Exception,e:
			print e
	return cldf.File(path,mode)
