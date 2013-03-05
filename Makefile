#WORK IN PRODRESS USE AT YOUR OWN RISK. NOT TESTED YET ON OTHER INFRASTRUCTURE THAN MACOS YOU HAVE BEEN WARNED !

# here are the thing that you should modify 

################################################################################################

# set your prefix to where you want to install clik.
# default is to let it in the current directory
PREFIX := $(shell pwd)

# set the path of the cfitsio lib. 
CFITSIOPATH := /usr/local
#CFITSIOPATH := /softs/cfitsio/3.24

#define your compilers and stuff
CC = gcc
FC = ifort

# ifort
# if you are using ifort set here where its lib are installed
# and check the runtime libs

# on my mac I got
IFORTLIBPATH = /usr/bin/ifort-2011-base/compiler/lib
IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifportmt -lifcoremt -lpthread

# on a linux machine, ifort 11.1
#IFORTLIBPATH = /softs/intel/fce/11.1.075/lib/intel64
#IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifport -lifcoremt -lpthread

# gfortran
# if you are using gfortran set here where the lib are installed
# and check the runtime libs
GFORTRANLIBPATH = /usr/lib
GFORTRANRUNTIME = -L$(GFORTRANLIBPATH) -lgfortran -lgomp

# if you are on linux and using mkl, you need to set this 
MKLROOT = /softs/intel/mkl/10.2.6.038/
# on mkl 10.3
#LAPACKLIBPATHMKL = -L$(MKLROOT)/lib/intel64
# on mkl 10.2
LAPACKLIBPATHMKL = -L$(MKLROOT)/lib/em64t

# pretty colors (comment to remove pretty colors)
COLORS = 1

# set the variable to the python cli to compile and install the python tools
PYTHON = python

################################################################################################

# you should not need to modify anything below


#temporary dirs
BDIR := $(shell pwd)/buildir
ODIR := $(shell pwd)/buildir/tmp

# tools
LD = gcc
INSTALL = install
ECHO = echo

# get the os
UNAME := $(shell uname -s)

ifeq ($(UNAME),Darwin)
OS = macos
else
OS = linux
endif

#defines for macos
SOMACOS = dylib
LIBPATHNAMEMACOS = DYLD_LIBRARY_PATH
#defines for linux
SOLINUX = so
LIBPATHNAMELINUX = LD_LIBRARY_PATH

ifeq ($(OS),macos)
SO = $(SOMACOS)
LIBPATHNAME = $(LIBPATHNAMEMACOS)
else
SO = $(SOLINUX)
LIBPATHNAME = $(LIBPATHNAMELINUX)
endif

#ifort
IFORTMODULEPATH = -module

#gfortran
GFORTRANMODULEPATH = -J

# this picks either ifort or gfortran, change those lines to set FRUNTIME and FMODULEPATH for your special case
ifeq ($(FC),ifort)
FLIBPATH = $(IFORTLIBPATH)
FRUNTIME = $(IFORTRUNTIME)
FMODULEPATH = $(IFORTMODULEPATH)
else
FLIBPATH = $(GFORTRANLIBPATH)
FRUNTIME = $(GFORTRANRUNTIME)
FMODULEPATH = $(GFORTRANMODULEPATH)
endif


# some defines (shared, relocatable openmp, etc)
CFPIC = -fPIC
COPENMP = -fopenmp

FFPIC = -fPIC
FOPENMP = -openmp

# check here that the SHARED variable contain the correct invocation for your CC
ifeq ($(OS),macos)
SHARED = -dynamiclib
else
SHARED = -shared -Bdynamic
endif

# get version of the code from the svn version
VERSION = $(strip $(shell cat svnversion)) MAKEFILE

# some more defines
#macos
DEFINESMACOS = -D HAS_RTLD_DEFAULT
#linux
DEFINESLINUX = 

DEFINESCOMMON = -D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D 'CLIKSVNVERSION="$(VERSION)"'


ifeq ($(OS),macos)
DEFINES = $(DEFINESMACOS) $(DEFINESCOMMON)
CM64 = -arch x86_64
FM64 = -arch x86_64
else
DEFINES = $(DEFINESLINUX) $(DEFINESCOMMON)
CM64 = -m64
FM64 = -m64
endif

INCLUDES = -I$(CFITSIOPATH)/include

# final CFLAG and FFLAGS
CFLAGS = $(CM64) $(COPENMP) $(CFPIC) $(DEFINES) -I src -I src/cldf -I src/minipmc -I src/lenslike/plenslike $(INCLUDES)
FFLAGS = $(FM64) $(FOPENMP) $(FFPIC) $(DEFINES) $(FMODULEPATH) $(ODIR)


# Lapack section

#macos I advise you to use the builtin blas lapack that are reasonnably efficient
LAPACKLIBPATHMACOS = /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current
LAPACKMACOS = -L$(LAPACKLIBPATHMACOS) -lBLAS -lLAPACK

# mkl I am assuming that the env variable MKLROOT contains the MKL root path
# if not define it here

LAPACKMKLCORELIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LAPACKMKL = -L$(LAPACKLIBPATHMKL) $(LAPACKMKLCORELIB)  -liomp5 -lpthread -lm

LAPACK_FUNC := dtrsv  dpotrf  dpotrs  dpotri  dtrtri  dtrmm  dtrmv  dgeqrf  dormqr  dsyev  dgesvd  dsymv  dgemv  dgemm  dsyrk  dsyr2k  daxpy  dtrsm  dsymm  dsyr  ddot
MKL_TO_INCLUDE := $(addprefix -u ,$(addsuffix _,$(LAPACK_FUNC)))
MKL_LIB_FULLPATH := $(filter $(addsuffix .a,$(addprefix %/lib,$(subst -l,,$(filter -l%,$(LAPACKMKLCORELIB))))),$(wildcard $(subst -L,,$(filter -L%,$(LAPACKMKL)))/lib*.a))

# pick lapack version
ifeq ($(OS),macos)
#macos lapack
LAPACK = $(LAPACKMACOS)
LAPACKLIBPATH = $(LAPACKLIBPATHMACOS)
LAPACKDEP =
else
#mkl !
LAPACK = -L$(BDIR) -llapack_clik
LAPACKLIBPATH = $(LAPACKLIBPATHMKL)
LAPACKDEP = $(BDIR)/liblapack_clik.$(SO)
endif

#if you want to point to your own version of lapack set the following variables
#LAPACK = -L/some/path -lsomefortranlapack -lsomedependencyforyourlapack
#LAPACKLIBPATH = /some/path
# leave this one empty
#LAPACKDEP = 

#CFITSIO link
CFITSIO =  -L$(CFITSIOPATH)/lib -lcfitsio

#final LDFLAG
LDFLAG = $(CM64) $(CFITSIO) $(LAPACK) $(FRUNTIME) -ldl -lm -lpthread

# define some path to find the codes
vpath %.c %f90 %.F90 src src/minipmc src/cldf src/CAMspec src/component_plugin/basic src/lenslike/plenslike
vpath %f90  src src/minipmc src/cldf src/CAMspec src/gibbs src/act_spt src/lowlike
vpath  %.F90 src src/minipmc src/cldf src/CAMspec src/gibbs src/act_spt src/lowlike

# define color output if needed
ifeq ($(COLORS),1)
NO_COLOR=\x1b[0m
GREEN_COLOR=\x1b[32;11m
RED_COLOR=\x1b[31;01m
BLUE_COLOR=\x1b[35;11m
endif

# all the code
TOOLS := $(addprefix $(ODIR)/,errorlist.o io.o distribution.o cldf.o)
CLIKMAIN := $(addprefix $(ODIR)/,clik.o lklbs.o lowly_common.o clik_helper.o)
CLIKLKL := $(addprefix $(ODIR)/,clik_lowlike.o clik_actspt.o clik_gibbs.o clik_CAMspec.o)
LENSLKL := $(addprefix $(ODIR)/,plenslike_dat_mono.o plenslike_dat_quad.o qest.o wignerd.o)
ACTSPTLKL := $(addprefix $(ODIR)/,Highell_options.f90.o Highell_subroutines.f90.o  Foregrounds_loading.f90.o ACT_equa_likelihood.f90.o SPT_reichardt_likelihood.f90.o ACT_south_likelihood.f90.o  SPT_keisler_likelihood.f90.o  Highell_likelihood.f90.o clik_actspt.f90.o)
CAMSPECLKL := $(addprefix $(ODIR)/,CAMspec.f90.o clik_CAMspec.f90.o)
LOWLIKELKL := $(addprefix $(ODIR)/,healpix_types.f90.o read_archive_map.f90.o read_fits.f90.o br_mod_dist.f90.o Planck_options.f90.o  Planck_teeebb_pixlike.f90.o  Planck_likelihood.f90.o clik_lowlike.f90.o)
GIBBSLKL := $(addprefix $(ODIR)/,comm_br_mod.f90.o clik_gibbs.f90.o)
CLIKLKL_F90:= $(ACTSPTLKL) $(CAMSPECLKL) $(GIBBSLKL) $(LOWLIKELKL)

CLIKLIB := $(TOOLS) $(CLIKMAIN) $(CLIKLKL) $(CLIKLKL_F90) $(LENSLKL) $(LAPACKDEP)


all: $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(BDIR)/clik_example_C $(BDIR)/clik_example_f90

install_dir: 
	@mkdir -p $(PREFIX)/bin
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@mkdir -p $(PREFIX)/share/clik

install: $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(BDIR)/clik_example_C $(BDIR)/clik_example_f90 $(LAPACKDEP) $(BDIR)/clik_profile.sh $(BDIR)/clik_profile.csh | install_dir
	@$(ECHO) "install libs $(BLUE_COLOR)libclik.$(SO) libclik_f90.$(SO)$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/lib $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(LAPACKDEP) $(PREFIX)/lib
	@$(ECHO) "install includes $(BLUE_COLOR)clik.h clik.mod$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/include $(NO_COLOR)"
	@$(INSTALL)  src/clik.h src/minipmc/maths_base.h src/minipmc/errorlist.h src/minipmc/io.h src/lapack_clik.h src/minipmc/pmc.h $(ODIR)/clik.mod $(PREFIX)/include
	@$(ECHO) "install clik_profile $(BLUE_COLOR)clik_profile.sh clik_profile.csh$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/clik_profile.sh $(BDIR)/clik_profile.csh $(PREFIX)/bin
	@$(ECHO) "install exec tools $(BLUE_COLOR)clik_example_C clik_example_f90$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/clik_example_C $(BDIR)/clik_example_f90 $(PREFIX)/bin


ifdef PYTHON
PYTHONPATH = $(PREFIX)/lib/`$(PYTHON) -c"import sys;print 'python%s/site-packages'%sys.version[0:3]"`
PYTHONEXE := `which $(PYTHON)`
else
PYTHONPATH := 
endif

$(BDIR)/clik_profile.sh: src/clik_profile.sh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s@CFITSIOLIBPATH@$(CFITSIOPATH)/lib@g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g;s!MPYTHONPATH!$(PYTHONPATH)!g" <$< >$@

$(BDIR)/clik_profile.csh: src/clik_profile.csh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s@CFITSIOLIBPATH@$(CFITSIOPATH)/lib@g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g;s!MPYTHONPATH!$(PYTHONPATH)!g" <$< >$@

	
$(BDIR):
	@mkdir $(BDIR)

$(ODIR): | $(BDIR)
	@mkdir $(ODIR)

$(CLIKLIB): | $(ODIR) $(ODIR)/.print_info

$(BDIR)/libclik.$(SO): $(CLIKLIB) 
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD)  $(SHARED)  $(LAPACK) $(LDFLAG) $^ -o $@

$(BDIR)/libclik_f90.$(SO): $(BDIR)/libclik.$(SO) $(addprefix $(ODIR)/,clik_fortran.o clik.f90.o)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD) $(SHARED)  $(LDFLAG) $(LAPACK) -L$(BDIR) -lclik $^ -o $@

$(BDIR)/clik_example_C: $(ODIR)/clik_example_c.o $(BDIR)/libclik.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(CC) $(LDFLAG) $(LAPACK) -L$(BDIR) -lclik $< -o $@

$(BDIR)/clik_example_f90: $(ODIR)/clik_example_f90.f90.o $(BDIR)/libclik_f90.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(FC) $(LDFLAG) $(LAPACK)  -L$(BDIR) -lclik_f90 -lclik $< -o $@

$(BDIR)/liblapack_clik.$(SO): |$(BDIR)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR),\n(see chapter 5 in http://software.intel.com/sites/products/documentation/hpc/mkl/lin/)\nusing the following command line:"
	gcc $(SHARED)  $(MKL_TO_INCLUDE) -Wl,--start-group $(MKL_LIB_FULLPATH) -Wl,--end-group -L$(IFORTLIBPATH) -L/lib -L/lib64 -liomp5 -lpthread -lm -o $@
	
$(ODIR)/%.o : %.c 
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(CC) -c $(CFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.f90 
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.F90 
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)

$(ODIR)/%.py: src/python/%.py
	@sed "s@PYTHONEXE@$(PYTHONEXE)@g;s@REPLACEPATH@$(PYTHONPATH)@g" <$< >$@
	@$(INSTALL) $@ $(PREFIX)/bin/$(subst .py,,$(@F))

$(ODIR)/.print_info: |$(ODIR)
	@$(ECHO) "\n$(BLUE_COLOR)Compile$(NO_COLOR) clik $(VERSION) "
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) CC = $(CC)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) FC = $(FC)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) CFLAGS = $(CFLAGS)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) FFLAGS = $(FFLAGS)"
	@$(ECHO) "$(BLUE_COLOR)Using the following lapack link line:$(NO_COLOR) $(LAPACK)"
	@$(ECHO) "$(BLUE_COLOR)Using the following cfitsio link line:$(NO_COLOR) $(CFITSIO)"
	@$(ECHO) "$(BLUE_COLOR)Using the following fortran runtime link line:$(NO_COLOR) $(FRUNTIME)"
	@$(ECHO) "$(BLUE_COLOR)Build dir:$(NO_COLOR) $(BDIR)"
	@$(ECHO)
	@touch $(ODIR)/.print_info

install_python: install $(addprefix $(ODIR)/, clik_add_free_calib.py clik_explore_1d.py prepare_actspt.py clik_get_selfcheck.py clik_example_py.py clik_join.py clik_disjoin.py clik_print.py prepare_wmap.py clik_extract_external.py) |$(ODIR)
	@LINK_CLIK="$(LDFLAG) $(LAPACK) -L$(PREFIX)/lib -lclik " $(PYTHON) setup.py build --build-base=$(ODIR) install --install-lib=$(PYTHONPATH)





clean:
	@$(ECHO) "$(BLUE_COLOR)Removing all in $(BDIR)$(NO_COLOR)"
	@rm -rf $(BDIR)

.PHONY :clean  LAPACK_PRINT LAPACK_DEP
