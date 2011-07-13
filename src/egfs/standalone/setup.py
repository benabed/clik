from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os.path as osp

import numpy as nm
import os

# remplace la chaine de charactere par le resultat de 
# gfortran -print-file-name=libgfortran.so 
# (ou sur mac
# gfortran -print-file-name=libgfortran.dylib )

fruntime = "/usr/local/gfortran/lib/gcc/i386-apple-darwin8.10.1/4.5.0/../../../libgfortran.dylib"


Lfruntime = osp.split(fruntime)[0]
Lfruntime = osp.realpath(Lfruntime)
lfruntime = "gfortran"

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("egfs", ["egfs.pyx","clik_egfs.c","distribution.c","io.c","errorlist.c"],
                             include_dirs = [nm.get_include()],
                             library_dirs=['.',Lfruntime],
                             libraries=["egfs",lfruntime]
                             )],
)