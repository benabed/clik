# set your prefix to where you want to install clik.
# default is to let it in the current directory
PREFIX := $(shell pwd)

#temporary dirs
BDIR := $(shell pwd)/buildir
ODIR := $(shell pwd)/buildir/tmp

#define your compilers and stuff
CC = gcc
FC = ifort
LD = gcc
INSTALL = install
ECHO = echo

UNAME := $(shell uname -s)

ifeq ($(UNAME),Darwin)
OS = macos
else
OS = linux
endif

#macos
SOMACOS = dylib
LIBPATHNAMEMACOS = DYLD_LIBRARY_PATH
#linux
SOLINUS = so
LIBPATHNAMELINUX = LD_LIBRARY_PATH

ifeq ($(OS),macos)
SO = $(SOMACOS)
LIBPATHNAME = $(LIBPATHNAMEMACOS)
else
SO = $(SOLINUX)
LIBPATHNAME = $(LIBPATHNAMELINUX)
endif

#ifort
# if you are using ifort set here where its lib are installed
IFORTLIBPATH = /usr/bin/ifort-2011-base/compiler/lib
# and check the runtime libs
IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifportmt -lifcoremt -lpthread
IFORTMODULEPATH = -module

#gfortran
#if you are using gfortran set here where the lib are installed
GFORTRANLIBPATH = /usr/lib
# and check the runtime libs
GFORTRANRUNTIME = -L$(GFORTRANLIBPATH) -lgfortran -lgomp
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


CFPIC = -fPIC
COPENMP = -fopenmp
CM64 = -arch x86_64

FFPIC = -fPIC
FOPENMP = -openmp
FM64 = -arch x86_64

VERSION = $(strip $(shell cat svnversion)) MAKEFILE
#macos
DEFINESMACOS = -D HAS_RTLD_DEFAULT
#linux
DEFINESLINUX = 

DEFINESCOMMON = -D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D 'CLIKSVNVERSION="$(VERSION)"'


ifeq ($(OS),macos)
DEFINES = $(DEFINESMACOS) $(DEFINESCOMMON)
else
DEFINES = $(DEFINESLINUX) $(DEFINESCOMMON)
endif


CFLAGS = $(CM64) $(COPENMP) $(CFPIC) $(DEFINES) -I src -I src/cldf -I src/minipmc -I src/lenslike/plenslike 
FFLAGS = $(FM64) $(FOPENMP) $(FFPIC) $(DEFINES) $(FMODULEPATH) $(ODIR)

# check here that the SHARED variable contain the correct invocation for your CC
ifeq ($(OS),macos)
SHARED = -dynamiclib
else
SHARED = -shared -Bdynamic
endif


# Lapack section

#macos I advise you to use the builtin blas lapack that are reasonnably efficient
LAPACKLIBPATHMACOS = /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current
LAPACKMACOS = -L$(LAPACKLIBPATHMACOS) -lBLAS -lLAPACK

# mkl I am assuming that the env variable MKLROOT contains the MKL root path
# if not define it here
#MKLROOT = /opt/intel/mkl

#mkl 10.3
LAPACKLIBPATHMKL103 = -L$(MKLROOT)/lib/intel64
LAPACKMKL103 = -L$(LAPACKLIBPATHMKL103) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core   -liomp5 -lpthread -lm
#mkl 10.2
LAPACKLIBPATHMKL102 = -L$(MKLROOT)/lib/emt64
LAPACKMKL102 = -L$(LAPACKLIBPATHMKL102) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core   -liomp5 -lpthread -lm
#mkl 10.1
LAPACKLIBPATHMKL101 = -L$(MKLROOT)/lib/emt64
LAPACKMKL101 = -L$(LAPACKLIBPATHMKL101) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core   -liomp5 -lpthread -lm

#set LAPACKMKL to the correct link incantation, here I am assuming that you are using LAPACK 10.3
LAPACKMKL = $(LAPACKMKL103)
LAPACKLIBPATHMKL = $(LAPACKLIBPATHMKL103)


LAPACK_FUNC := dtrsv  dpotrf  dpotrs  dpotri  dtrtri  dtrmm  dtrmv  dgeqrf  dormqr  dsyev  dgesvd  dsymv  dgemv  dgemm  dsyrk  dsyr2k  daxpy  dtrsm  dsymm  dsyr  ddot
MKL_TO_INCLUDE := $(addprefix -u ,$(addsuffix _,$(LAPACK_FUNC)))
MKL_LIB_FULLPATH := $(filter $(addsuffix .a,$(addprefix %/lib,$(subst -l,,$(filter -l%,$(LAPACKMKL))))),$(wildcard $(subst -L,,$(filter -L%,$(LAPACKMKL)))/lib*.a))


ifeq ($(OS),macos)
#macos lapack
LAPACK = $(LAPACKMACOS)
LAPACKLIBPATH = $(LAPACKLIBPATHMACOS)
LAPACKDEP =
else
#mkl !
LAPACK = $(LAPACKMKL) 
LAPACKLIBPATH = $(LAPACKLIBPATHMKL)
LAPACKDEP = $(BIDR)/lapack_clik.$(SO)
endif

CFITSIO = -L/usr/local/lib -lcfitsio

LDFLAG = $(CM64) $(CFITSIO) $(LAPACK) $(FRUNTIME) -ldl -lm -lpthread

vpath %.c %f90 %.F90 src src/minipmc src/cldf src/CAMspec src/component_plugin/basic src/lenslike/plenslike
vpath %f90  src src/minipmc src/cldf src/CAMspec src/gibbs src/act_spt src/lowlike
vpath  %.F90 src src/minipmc src/cldf src/CAMspec src/gibbs src/act_spt src/lowlike


NO_COLOR=\x1b[0m
GREEN_COLOR=\x1b[32;01m
RED_COLOR=\x1b[31;01m
BLUE_COLOR=\x1b[32;11m


TOOLS := $(addprefix $(ODIR)/,errorlist.o io.o distribution.o cldf.o)
CLIKMAIN := $(addprefix $(ODIR)/,clik.o lklbs.o lowly_common.o smica.o clik_helper.o)
CLIKLKL := $(addprefix $(ODIR)/,clik_lowlike.o clik_actspt.o clik_gibbs.o clik_CAMspec.o clik_hfipack.o clik_parametric.o clik_parametric_addon.o basic.o)
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
	$(INSTALL) -C $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(LAPACKDEP) $(PREFIX)/lib
	@$(ECHO) "install includes $(BLUE_COLOR)clik.h clik.mod$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/include $(NO_COLOR)"
	$(INSTALL) -C src/clik.h src/minipmc/errorlist.h src/minipmc/io.h src/lapack_clik.h src/minipmc/pmc.h $(ODIR)/clik.mod $(PREFIX)/include
	@$(ECHO) "install clik_profile $(BLUE_COLOR)clik_profile.sh clik_profile.csh$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	$(INSTALL) -C $(BDIR)/clik_profile.sh $(BDIR)/clik_profile.csh $(PREFIX)/bin
	@$(ECHO) "install exec tools $(BLUE_COLOR)clik_example_C clik_example_f90$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	$(INSTALL) -C $(BDIR)/clik_example_C $(BDIR)/clik_example_f90 $(PREFIX)/bin


$(BDIR)/clik_profile.sh: src/clik_profile.sh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g" <$< >$@
$(BDIR)/clik_profile.csh: src/clik_profile.csh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g" <$< >$@
	
$(BDIR):
	mkdir $(BDIR)

$(ODIR): | $(BDIR)
	mkdir $(ODIR)

$(CLIKLIB): | $(ODIR)

$(BDIR)/libclik.$(SO): $(CLIKLIB)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD)  $(SHARED)  $(LDFLAG) $^ -o $@

$(BDIR)/libclik_f90.$(SO): $(BDIR)/libclik.$(SO) $(addprefix $(ODIR)/,clik_fortran.o clik.f90.o)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD) $(SHARED)  $(LDFLAG) -L$(BDIR) -lclik $^ -o $@

$(BDIR)/clik_example_C: $(ODIR)/clik_example_c.o $(BDIR)/libclik.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(CC) $(LDFLAG) -L$(BDIR) -lclik $< -o $@

$(BDIR)/clik_example_f90: $(ODIR)/clik_example_f90.f90.o $(BDIR)/libclik_f90.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(FC) $(LDFLAG) -L$(BDIR) -lclik_f90 -lclik $< -o $@

$(BIDR)/lapack_clik.$(SO): |$(BDIR)
	gcc $(shared)  $(MKL_TO_INCLUDE) -Wl,--start-group $(MKL_LIB_FULLPATH) -Wl,--end-group $(FRUNTIME) -L/lib -L/lib64 -liomp5 -lpthread -lm -o $@

$(ODIR)/%.o : %.c 
	@$(ECHO) "$(BLUE_COLOR)$< $(NO_COLOR) -> $(BLUE_COLOR) $(@) $(NO_COLOR)"
	@$(CC) -c $(CFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.f90 
	@$(ECHO) "$(BLUE_COLOR)$< $(NO_COLOR) -> $(BLUE_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.F90 
	@$(ECHO) "$(BLUE_COLOR)$< $(NO_COLOR) -> $(BLUE_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)


clean:
	rm -rf $(BDIR)

.PHONY :clean