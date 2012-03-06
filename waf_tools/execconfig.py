import waflib.TaskGen
import waflib.Task as Task
from waflib import Utils
import waflib
import os.path as osp
import os

def uniqify(lst):
  rlst = lst
  for v in lst:
    if v in rlst:
      continue
    rlst.append(v)
  return rlst
def ptrquote(st):
  res = ""
  for v in st:
    if v=='"':
      res +='\\"'
    else:
      res+=v
  return res

@waflib.TaskGen.feature("build_pkgconfig")
def build_pkgconfig(self):
  from waflib.Tools.ccroot import USELIB_VARS
  if self.flavor=='c':  
    USELIB_VARS['build_pkgconfig']   = set(['INCLUDES', 'DEFINES', 'CPPFLAGS', 'CFLAGS']+['LIB', 'STLIB', 'LIBPATH', 'STLIBPATH', 'LINKFLAGS', 'RPATH', 'LINKDEPS'])
    cf = ['CPPFLAGS', 'CFLAGS']
  else:
    USELIB_VARS['build_pkgconfig']   =set(['FCFLAGS','DEFINES','INCLUDES']+['LIB','STLIB','LIBPATH','STLIBPATH','LINKFLAGS','RPATH','LINKDEPS'])
    cf = ['FCFLAGS']

  #USELIB_VARS['cprogram']
  self.process_use()
  self.propagate_uselib_vars()
  vrs = dict([(v,list(set(self.env[v]))) for v in USELIB_VARS['build_pkgconfig']])

  includepath = ptrquote(" ".join([self.env.CPPPATH_ST%v for v in uniqify(vrs["INCLUDES"])]))
  libpath = ptrquote(" ".join([self.env.LIBPATH_ST%v for v in uniqify(vrs["LIBPATH"])]))
  rlibpath = ptrquote(" ".join([self.env.RPATH_ST%v for v in uniqify(vrs["RPATH"])]))
  stlibpath = ptrquote(" ".join([self.env.LIBPATH_ST%v for v in uniqify(vrs["STLIBPATH"])]))
  libs = ptrquote(" ".join([self.env.LIB_ST%v for v in uniqify(vrs["LIB"])]))
  stlibs = ptrquote(" ".join([self.env.STLIB_ST%v for v in uniqify(vrs["STLIB"])]))
  defines = ptrquote(" ".join([self.env.DEFINES_ST%v for v in uniqify(vrs["DEFINES"])]))
  cfs = []
  for tt in cf+["LINKFLAGS"]:
    cfs += vrs[tt]
  cflags = ptrquote(" ".join(uniqify(cfs)))

  #print includepath
  #print libpath
  #print rlibpath
  #print stlibpath
  #print libs
  #print stlibs
  #print cflags
  #print defines
  
  alibs = ""
  if libs:
    alibs += (self.env.SHLIB_MARKER or "") +" ".join([rlibpath,libpath,libs])
  if stlibs:
    alibs += (self.env.STLIB_MARKER or "") +" ".join([srlibpath,stlibs])

  f=open(osp.join(self.env.BINDIR,self.target),"w")
  print >>f,config_tpl%(" ".join((includepath,defines,cflags)),alibs)
  f.close()  
  os.chmod(osp.join(self.env.BINDIR,self.target),Utils.O755)

config_tpl = """#! /usr/bin/env python
# don't do much for now
from optparse import OptionParser
parser = OptionParser()

parser.add_option("--cflags", action="store_true",
                  help="only the cflags")
parser.add_option("--libs", action="store_true",
                  help="only libflags")

(options, args) = parser.parse_args()


res={}
cflags = "%s"
libs = "%s"

if (not options.libs) and (not options.cflags):
  options.libs=True
  options.cflags=True

if options.cflags:
  print cflags,
if options.libs:
  print libs,
print

"""