#try to support many flavours of lapack
import autoinstall_lib as atl
from waflib import Logs
import os.path as osp

version = "lapack-3.3.1"
tool = "lapack-3.3.1"
lapack_funcs = "dtrsv dpotrf dpotri dtrtri dtrmm dtrmv dgeqrf dormqr dsyev dgesvd dsymv dgemv dgemm dsyrk dsyr2k daxpy dtrsm dsymm dsyr ddot"     

def options(ctx):
  atl.add_lib_option("lapack",ctx,install=True)
  grp = ctx.parser.get_option_group("--lapack_install")
  grp.add_option("--lapack_mkl",action="store",default="",help="if lapack is mkl, location of the mkl install")

def configure(ctx):
  # no veclib !
  if "vecLib" in ctx.options.lapack_lib or "veclib" in ctx.options.lapack_include:
    #must be darwin !
    if ctx.env.has_lapack:
      raise Exception("Accelerate framework is not supported at this time")
    else:
      Logs.pprint("RED","Apple vecLib not supported")
      return
  
  lapack_extradefs = ["HAS_LAPACK"]
  lapack_libs = ["BLAS","LAPACK"]
  lapack_includes = ["lapack.h","blas.h"]

  if "mkl" in ctx.options.lapack_lib.lower() or "mkl" in ctx.options.lapack_include.lower() or "mkl" in ctx.options.lapack_link or ctx.options.lapack_mkl:
    ctx.env.mkl = True
    lapack_extradefs += ["HAS_MKL"]
    lapack_includes = ["mkl_lapack.h","mkl_blas.h"]
    if ctx.options.lapack_mkl:
      if ctx.env.has_ifort==False:
        raise Exception("cannot use MKL without ifort")
      if "framework" in ctx.options.lapack_mkl.lower():
        # guess we are on macosx
        # get the path of the framework
        if ctx.options.lapack_mkl[-1] == "/":
          fpath,fname = osp.split(ctx.options.lapack_mkl[:-1])
        else:
          fpath,fname = osp.split(ctx.options.lapack_mkl)
        fname = fname.split(".")[0]
        ctx.options.lapack_include =  ctx.options.lapack_mkl+"/Headers"
        ctx.options.lapack_lib =  ctx.options.lapack_mkl+"/Libraries/universal"
        if ctx.options.lapack_link=="":
          ctx.options.lapack_link = "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
      else:
        # assume it's 10 on linux
        # check whether it's 10.3
        if ctx.options.m32:
          libsuffix="/lib/32"
          libdep = "-lmkl_intel"
        else:
          libsuffix="/lib/em64t"
          libdep = "-lmkl_intel_lp64"
        if ctx.options.lapack_link=="":
          ctx.options.lapack_link = "-lmkl_lapack -lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread -lmkl_def" + libdep
        if not ctx.options.m32 and osp.exists(ctx.options.lapack_mkl+"/lib/intel64"):
          libsuffix="/lib/intel64"
          ctx.options.lapack_link = "-lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread -lmkl_def" + libdep
        ctx.options.lapack_include=ctx.options.lapack_mkl+"/include"
        ctx.options.lapack_lib=ctx.options.lapack_mkl+libsuffix+":".join([""]+ctx.env.LIBPATH_fc_runtime)
  iall = atl.shouldIinstall_all(ctx,"lapack")
  if ctx.options.lapack_install or ctx.options.lapack_islocal or ctx.options.lapack_forceinstall or iall:
    ctx.env.append_value("LIBPATH_lapack",ctx.env.LIBPATH_fc_runtime)
    ctx.env.append_value("RPATH_lapack",ctx.env.RPATH_fc_runtime)
    ctx.env.append_value("LIB_lapack",ctx.env.LIB_fc_runtime)
    lapack_libs = ["lapack_clik","blas_clik"]
    lapack_includes = ["lapack_clik.h"]
    lapack_extradefs += ["LAPACK_CLIK"]
    
  atl.conf_lib(ctx,"lapack",lapack_libs,lapack_funcs.split(),lapack_includes,defines=lapack_extradefs,install=installlapack)

def installlapack(ctx):
  filen = version+".tgz"
  atl.installsmthg_pre(ctx,"http://www.netlib.org/lapack/"+filen,filen)
  from waflib import Utils,Errors
  dii = {"FCC":ctx.env.FC,"FCFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fcshlib),"FLINKFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.LINKFLAGS_fcshlib),"SO":ctx.env.shsuffix,"MFLAG":" ".join(ctx.env.FCFLAGS) }
  Logs.pprint("PINK","build blas")
  f=open("build/%s/make.inc"%version,"w")
  print >>f,make_inc_blas%dii
  f.close()
  cmdline = "cd build/%s; make -j blaslib"%version
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%version)
  Logs.pprint("PINK","build lapack")
  f=open("build/%s/make.inc"%version,"w")
  print >>f,make_inc_lapack%dii
  f.close()
  cmdline = "cd build/%s; make -j lapacklib"%version
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%version)
  
  import shutil
  shutil.copyfile("build/%s/liblapack_clik.%s"%(version,ctx.env.shsuffix), osp.join(ctx.env.LIBDIR,"liblapack_clik.%s"%ctx.env.shsuffix))
  shutil.copyfile("build/%s/libblas_clik.%s"%(version,ctx.env.shsuffix), osp.join(ctx.env.LIBDIR,"libblas_clik.%s"%ctx.env.shsuffix))

  f=open(osp.join(ctx.env.PREFIX,"include/lapack_clik.h"),"w")
  for fnc in lapack_funcs.split():
    print >>f,"#define %s %s_"%(fnc,fnc)
  print >>f,extra_inc
  f.close()
  
make_inc_lapack="""
SHELL = /bin/sh
FORTRAN  = %(FCC)s %(FCFLAGS)s
OPTS     =
DRVOPTS  = $(OPTS)
NOOPT    = -g -O0
TIMER    = INT_CPU_TIME
LOADER   = %(FCC)s
LOADOPTS = %(MFLAG)s

BLASLIB      = ../../libblas_clik.%(SO)s
ARCH = %(FCC)s 
ARCHFLAGS = %(FLINKFLAGS)s -L../ -lblas_clik -o
RANLIB = echo
LAPACKLIB    = liblapack_clik.%(SO)s
"""

make_inc_blas="""
SHELL = /bin/sh
FORTRAN  = %(FCC)s %(FCFLAGS)s
OPTS     =
DRVOPTS  = $(OPTS)
NOOPT    = -g -O0
TIMER    = INT_CPU_TIME

BLASLIB      = ../../libblas_clik.%(SO)s
ARCH = %(FCC)s 
ARCHFLAGS = %(FLINKFLAGS)s -o
RANLIB = echo
LAPACKLIB    = liblapack_clik.%(SO)s
"""

extra_inc = """
void dtrsv(const char *uplo, const char *trans, const char *diag, const int  *n,
           const double *a, const int *lda, double *x, const int *incx);
void dpotrf( char* uplo, int * n, double* a, int * lda, int * info );
void dpotri( char* uplo, int * n, double* a, int * lda, int * info );
void dgemv(const char *trans, const int *m, const int *n, const double *alpha,
           const double *a, const int *lda, const double *x, const int *incx,
           const double *beta, double *y, const int *incy);
void dsyrk(const char *uplo, const char *trans, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *beta,
           double *c, const int *ldc);
void dsyr2k(const char *uplo, const char *trans, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
void dgesvd( char* jobu, char* jobvt, int * m, int * n, double* a, int * lda, double* s, double* u, int * ldu, double* vt, int * ldvt, double* work, int * lwork, int * info );
void dgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
void dtrtri( char* uplo, char* diag, int * n, double* a, int * lda, int * info );
void dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dtrmv(const char *uplo, const char *transa, const char *diag, const int *n,
           const double *a, const int *lda, double *b, const int *incx);
void dgeqrf( int * m, int * n, double* a, int * lda, double* tau, double* work, int * lwork, int * info );
void dormqr( char* side, char* trans, int * m, int * n, int * k, double* a, int * lda, double* tau, double* c, int * ldc, double* work, int * lwork, int * info );
void dsyev( char* jobz, char* uplo, int * n, double* a, int * lda, double* w, double* work, int * lwork, int * info );
void dsymv(const char *uplo, const int *n, const double *alpha, const double *a, const int *lda,
           const double *x, const int *incx, const double *beta, double *y, const int *incy);
void daxpy(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dsyr(const char *uplo, const int *n, const double *alpha, const double *x, const int *incx,
         double *a, const int *lda);
void dsymm(const char *side, const char *uplo, const int *m, const int *n,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
double ddot(int* N,double *DX, int* INCX,double *DY,int* INCY);
           
         
"""