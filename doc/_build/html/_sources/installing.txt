Installing
==========

The library has a lot of external dependencies. This should not evolve in the future, hopefully.

First, it relies heavily on functions that are provided in pmclib. So the first step is to compile pmclib. pmclib will require gsl (1.10+), fftw3, hdf5 (1.8+)  and the intel mkl. Apart from the mkl, all those can be installed for you by the pmclib installer.

Luckily most of the dependencies of clik are shared with pmclib, so after this step, there is only one extra (mandatory) dependency the c version of the Healpix library, compiled as a shared library.

Optionnaly, you will also need a fortran compiler for the f90 wrapper, and python, cython, numpy and h5py for the python wrapper and the synthetic smica likelihood generator.

I the following, 64 bits compilation is assumed. To switch to 32 bits, replace the --m64 command line option by --m32.

Installing pmclib
-----------------

Please only use the pmclib provided along with the clik lib. The installation in this package is simplified (sic...) compared to the one currently publicly available (1.1). That will change in the future.

Assuming you have unzipped and untard the pmclib package, and that the install of your mkl library is located at MYMKLPATH, the following line should configure your pmclib before building and install for you the gsl, fftw3 and hdf5 lib into at the path given in the prefix option (in this case, your current working directory). The library will also install itself in the directory pointed by prefix.

$ ./waf configure --m64 --prefix=`pwd` --lapack_mkl=MYMKLPATH --hdf5_install --gsl_install --fftw3_install

If you already have inatalled some of the dependency (say for example, gsl), replace the option --gsl_install by gsl_prefix=MYGSLPATH, so that MYGSLPATH/lib points to the libs, and MYGSLPATH/include points to the include files of your gsl install.

More options are available by running 

$ ./waf --help

Once the configure step is passed successfully, you can build and install the pmclib using the command

$ ./waf install

Installing clik
----------------

After the installation of pmclib, cliklib can also be installed. Assuming that you have unzipped and untard the clik package, 

$ ./waf configure --m64 --prefix=`pwd` --pmc_prefix=MYPMCLIBPATH --chealpix_prefix=MYHEALPIXPATH

If h5py cannot be found it can be installed adding the --h5py_install command line option to the previous configure command.

The library and its dependency will be installed at the path pointed by prefix (in this case the current directory).

Once the configure step is passed successfully, you can build and install the clik using the command

$ ./waf install



