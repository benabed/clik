Installing
==========

The library has a lot of external dependencies. 

Only a few packages or tools must be installed before installing clik. All the other can be downloaded and installed automatically by the clik installer.

Furthermore some of the dependency can be absent, and will only reduce the number of extra functionnalities provided by clik (like merging likelihood, or simulating them).

Requisites
----------

Mandatory requisites that cannot be installed automatically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Python <http://python.org>`_ (>2.5) to run the installer (waf), **modern c and fortran compilers** (icc, gcc, ifort anf gfortran are ok) are absolute requisites. 
Having a full python installation (i.e. including the library and header, often called *python-devel* or *python-dev* in package managers) is very strongly advised.

Mandatory requisites that will be installed automatically if absent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `gsl <http://www.gnu.org/software/gsl/>`_ (version>1.10), `hdf5 <http://www.hdfgroup.org/HDF5>`_ (version>1.8), `healpix c and f90 <http://healpix.jpl.nasa.gov/>`_ (version>2.20a shared) libraries with their dependency (i.e. cfitsio) as well as blas and lapack distribution (preferably intel MKL) are needed for the core functionalities of clik. If absent (or not available in the correct flavor), they can be installed automatically using the options described below.

Optional requisites 
^^^^^^^^^^^^^^^^^^^

A python distribution including the header and library (if absent, this will not be installed automatically), and the `h5py <http://alfven.org/wp/hdf5-for-python/>`_ (1.3<=version<2),  `numpy <http://numpy.scipy.org/>`_ (version>1.1) and `cython <http://cython.org/>`_ python package (version>1.12) are needed to provide the (optional) clik python wrapper and the simulation, merging , splitting and printing tools. The three python package will only be installed automatically if the python header and library are available on the system.

Special dependency : pmclib
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clik uses some of the facilities available in pmclib. Instead of requiring pmclib, those facilities are included in the clik distribution. This prevents a stock clik install to be linked with pmclib. An option is available to point the installer to a pmclib install, resolving this issue.


Install tool
------------

The library and executables must be installed with the 'waf' tool. It is distributed in the package. Please have a look at `the waf webpage <http://waf.googlecode.com>`_.

The installer must first be configured using::

    $> ./waf configure

then the codes can be compiled and installed with::

    $> ./waf install

The command line options given during the `configure` step allow to pass to the installer the locations of the different dependencies. For a complete list of options, do::

	$> ./waf configure --help

Only the most important ones will be described here.

Simplest install
^^^^^^^^^^^^^^^^

Setting the option ``--install_all_deps`` will install all the dependencies for you. This is the recommanded base option. Installing clik can be as simple as::

	$> ./waf configure --install_all_deps

Beware, this is slow and use a lot of disk space.

It can be desirable not to install all the dependencies if some are known to be available. Using any of the ``--XXX_prefix``, ``--XXX_islocal``, ``--XXX_installifneeded``, ``--XXX_include``, ``--XXX_lib``, ``XXX-link`` described below will cancel the installing option for the dependency ``XXX``.

The special dependency pmclib is not affected by this option.

Advanced install options
------------------------

Installing with a particular Python executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to install clik with a python install different from the default one. For example if the default python installation does not contains the required header and libraries. To do so, call waf this way::

    $> /path/to/special/python waf configure 

and then::

    $> /path/to/special/python waf install 


Bypassing the default compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To bypass the c compiler detection, set the ``CC`` environment variable. 
To bypass the fortran compiler detection, set the ``FC`` environment variable. Beware, you can only set the ``FC`` environment variable to either an intel fortran compiler or a gfortran compiler. 

Shortcuts for some classical cases are provided:

    * ``--icc`` causes the installer to use icc as c compiler.
    * ``--ifort`` causes the installer to use ifort as fortran compiler.
    * ``--gcc`` causes the installer to use gcc as c compiler.
    * ``--gfortran`` causes the installer to use gfortran as fortran compiler.


Setting the architecture
^^^^^^^^^^^^^^^^^^^^^^^^

The architecture (32 or 64bits) can be set using the ``--m32`` or ``--m64`` flags. 64bits is the default.

Setting installation path
^^^^^^^^^^^^^^^^^^^^^^^^^

The installation path can be set using the ``--prefix=SOMEPATH`` option. Default is to install in the current directory.


Telling clik to install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting the option ``--XXX_install`` will tell the clik install system to attempt at installing the dependency. It will first download it from the net, and then try to compile and install it. Installation will be perfomed at the same location as where clik itself will install. The download and install will only be done if needed. 
This is only possible for the following dependencies : gsl, hdf5, healpix, blas/lapack, numpy, cython, h5py. This option will install the code even if it is already available.
Using the option ``--XXX_install`` is recommanded.

To install only if a package is unavailable, use the option ``--XXX_installifneeded``. This last option is only recommanded to advanced users 


Setting the location of a library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Those options are only for advanced users.
The location of the library dependencies (gsl, hdf5, healpix, blas/lapack) must be known to the installer. By default, it will look for them in the classical system 
locations:  ``/usr/lib``, ``/usr/lib64``, ``/usr/local/lib``, ``/usr/local/lib64`` for the library, ``/usr/include`` and ``/usr/local/include`` for the include files. One can 
change the lookup path on a library by library basis. If a given dependency, ``XXX``, is installed on the system such that its lib are in ``SOMEPREFIXPATH/lib`` and its 
include files in ``SOMEPREFIXPATH/include``, setting the command line option ``--XXX_prefix=SOMEPREFIXPATH``  will allow the clik install system. If ``SOMEPREFIXPATH`` is identical to the the install path of clik, this option can be replaced by ``-XXX_islocal``.

If the library are at 
``SOMEWEIRDPATH`` and the includes at ``SOMEDIFFERENTPATH``, then setting the two options  ``--XXX_lib=SOMEWEIRDPATH --XXX_include=SOMEDIFFERENTPATH`` will allow the clik 
install system to find them.

Finally, if the name of the library files differs from the usual ones one can set the option ``--XXX_link=THELINKLINE``.

Using these options allow to point the installer to a pmclib install in order to allow the linking of clik with pmclib.


Special case: the mkl library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option is only for advanced users.
The blas/lapack distribution installed automatically is a very inefficient one. To improve the performance of clik (especially the low-l pixel based likelihood), one is advised to use the MKL library, which is fully supported and allow the use of shared memory computer architectures.
A special option is present to simplify the install using the intel MKL library: setting the option ``--lapack_mkl=PATH_OF_THE_MKL_INSTALL`` will allow clik to pick the correct set of libraries for the particular version of the mkl package (version 10.1, 10.2 and 10.3 tested).
Setting this option will cancel the ``--install_all_deps`` option for the lapack dependency only.

Special case: WMAP likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clik can provide a wrapper to the wmap7 likelihood. It need to now where the sources of the likelihood are located to compile against them. One must set the option ``--wmap_src=WMAP7SRCPATH`` or let the install system download it for you by setting the option ``--wmap_install``. Note that to actually use this likelihood, one must also download the data files and prepare clik likelihood files from them. Look at :ref:`WMAP`.


Special case: Healpix
^^^^^^^^^^^^^^^^^^^^^

Clik requires a specialy build healpix library. Namely, it insist on using a repositionnable (or better shared) version of the healpix library. 
This option is currently not available for the fortran version of the lib (as of version 2.20a). The clik installer know how to produce this special version
of healpix for you. Thus except if you really know what you are doing, and even if you already have healpix installed on your system, 
using the option ``--healpix_install`` is very strongly recommanded.

Putting it all together
^^^^^^^^^^^^^^^^^^^^^^^

The following command::

    $> ./waf configure --gsl_install --hdf5_install --lapack_install \
       --healpix_install --h5py_install --cython_install --numpy_install --wmap_install

will tell the clik install system to install all the possible external dependency in the current directory. All will be compiled in 64bit mode. Clik will be compiled in 64bit as well and installed in the same directory. This is the same as::

    $> ./waf configure --install_all_deps

The following command::

    $> ./waf configure --lapack_mkl=/opt/intel/mkl \
       --healpix_install --hdf5_install --h5py_install 

will tell the clik install system to install healpix, hdf5 and h5py. All the other dependency will be looked up in the classical locations. The blas/lapack library 
will be the one from an mkl install located at --lapack_mkl=/opt/intel/mkl. Clik will be compiled in 64bit and installed in the current directory.

 
Best advanced choice 
^^^^^^^^^^^^^^^^^^^^

Use a mkl lapack install and let the other dependencies on auto install::

    $> ./waf configure --lapack_mkl=/opt/intel/mkl --install_all_deps \
        --cython_installifneeded --numpy_installifneeded --gsl_installifneeded

This will use your mkl libraries from ``/opt/intel/mkl``, test if numpy, cython and gsl are installed on your computer (often the case) if not install them, 
and finally install all the other requirements (helpaix, hdf5 and its python wrapper).

Environment variables
---------------------

Depending of your shell, a configuration file named ``clik_profile.sh`` of ``clik_profile.csh`` will be installed in the ``bin`` directory at the install location of clik. One can source it on the command line, or include it in its startup configuration file to set the environment variable needed by clik.


