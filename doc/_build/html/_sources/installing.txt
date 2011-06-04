Installing
==========

The library has a lot of external dependencies. 

You will need `python <http://python.org>`_ (version>2.5, with library and header files) a c compiler (gcc or icc), a fortran compiler (gfortran or ifort), `gsl <http://www.gnu.org/software/gsl/>`_ (version>1.10), `hdf5 <http://www.hdfgroup.org/HDF5>`_ (version>1.8), `healpix c and f90 <http://healpix.jpl.nasa.gov/>`_ (version>2.20), a blas and lapack distribution (preferably intel MKL), the `h5py <http://alfven.org/wp/hdf5-for-python/>`_ python package (version>1.3), the `numpy <http://numpy.scipy.org/>`_ python package (version>1.1), the `cython <http://cython.org/>`_ python package (version>1.12). Optionnaly, it can use pmclib (version>1.5) if you have it on your system, and retrieve the location of most of those external library from it (as pmclib as mostly the same requirements).

If some of those dependencies are unavaliable, clik will lack a few facilities.
If the fortran compiler is absent, neither the BOPIX likelihood nor the clik fortran wrapper and example code will be compiled.
If the library and headr of the python package or any of the required modules are absent, neither the clik python wrapper, nor the likelihood file manipulation utilities will be compiled.

Finally, the clik installer can also (try to) install for you some of the dependencies, namelly, gsl, hdf5, blas/lapack, healpix, h5py, numpy and cython. In this case, no particular system specific optimization will be attempted. This is mainly a problem for the blas/lapack library that will be much less efficient than (for example) the intel MKL library, and will not take advantage of multicore cpus.

Install tool
------------

The library and executables must be installed with the 'waf' tool. It is distributed in the package. Please have a look at `the waf webpage <http://waf.googlecode.com>`_.

The installer must first be configured using::

    $> ./waf configure

then the codes can be compiled and installed with::

    $> ./waf install

The command line options given during the `configure` step allow to pass to the installer the locations of the different dependencies. For a complete list of options, do::

$> ./waf configure

Only the most important one will be described here.

Setting the architecture
^^^^^^^^^^^^^^^^^^^^^^^^

The architecture (32 or 64bits) can be set using the ``--m32`` or ``--m64`` flags. 64bits is the default.

Setting installation path
^^^^^^^^^^^^^^^^^^^^^^^^^

The installation path can be set using the ``--prefix=SOMEPATH`` option. To isntall in the current directory, one can use the ``--local`` option which is equivallent to ``--prefix=`pwd```


Setting the location of a library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The location of the library dependencies (gsl, hdf5, healpix, blas/lapack) must be known to the installer. By default, it will look for them in the classical system 
locations:  ``/usr/lib``, ``/usr/lib64``, ``/usr/local/lib``, ``/usr/local/lib64`` for the library, ``/usr/include`` and ``/usr/local/include`` for the include files. One can 
change the lookup path on a library by library basis. If a given dependency, ``XXX``, is installed on the system such that its lib are in ``SOMEPREFIXPATH/lib`` and its 
include files in ``SOMEPREFIXPATH/include``, setting the command line option ``--XXX_prefix=SOMEPREFIXPATH``  will allow the clik install system. If ``SOMEPREFIXPATH`` is identical to the the install path of clik, this option can be replaced by ``-XXX_islocal``.

If the library are at 
``SOMEWEIRDPATH`` and the includes at ``SOMEDIFFERENTPATH``, then setting the two options  ``--XXX_lib=SOMEWEIRDPATH --XXX_include=SOMEDIFFERENTPATH`` will allow the clik 
install system to find them.

Finally, if the name of the library files differs from the usual ones one can set the option ``--XXX_link=THELINKLINE``.


Special case: the mkl library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An extra special option is present for the blas/lapack dependency when one want to compile against the intel MKL library. Setting the option ``--lapack_mkl=PATH_OF_THE_MKL_INSTALL`` will allow clik to pick the correct set of libraries for the particular version of the mkl package (version 10.1, 10.2 and 10.3 tested).

Telling clik to install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting the option ``--XXX_install`` will tell the clik install system to attempt at installing the dependency. It will first download it from the net, and then try to compile and install it. Installation will be perfomed at the same location as where clik itself will install. The download and install will only be done if needed. 
This is only possible for the following dependencies : gsl, hdf5, healpix, blas/lapack, numpy, cython, h5py.


Putting it all together
^^^^^^^^^^^^^^^^^^^^^^^

The following command::

    $> ./waf configure --local --gsl_install --hdf5_install --lapack_install --healpix_install --h5py_install --cython_install --numpy_install

will tell the clik install system to install all the possible external dependency in the current directory. All will be compiled in 64bit mode. Clik will be compiled in 64bit as well and installed in the same directory.

The following command::

    $> ./waf configure --local --lapack_mkl=/opt/intel/mkl --healpix_install --hdf5_install --h5py_install 

will tell the clik install system to install healpix, hdf5 and h5py. All the other dependency will be looked up in the classical locations. The blas/lapack library 
will be the one from an mkl install located at --lapack_mkl=/opt/intel/mkl. Clik will be compiled in 64bit and installed in the current directory.

    
Environment variables
---------------------

Depending of your shell, a configuration file named ``clik_profile.sh`` of ``clik_profile.csh`` will be installed in the ``bin`` directory at the install location of clik. One can source it on the command line, or include it in its startup configuration file to set the environment variable needed by clik.


