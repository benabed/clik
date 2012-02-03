#include "fakedf.h"

int main(int argc, char **argv) {
  fdf *df;
  error *_err,**err;
  int a;
  double b;
  char *c;
  double *cr;
  int sz;
  int i;
  long *er;

  _err = NULL;
  err = &_err;
_DEBUGHERE_("",""); 
  df = fdf_open("test_fdf.clik",err);
  quitOnError(*err,__LINE__,stderr);
_DEBUGHERE_("%p",df);

  a = fdf_readint(df,"d/a",err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("a %d",a);

  b = fdf_readfloat(df,"b",err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("b %g",b);
  
  c = fdf_readstr(df,"c",err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("c '%s'",c);
  
  c = fdf_readstr(df,"d/b",err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("d/b '%s'",c);
  
  sz = 0;
  cr = fdf_readfloatarray(df,"d/c",&sz,err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("d/c size %d",sz);
  for(i=0;i<sz;i++) {
    _DEBUGHERE_("d/c[%d] %g",i,cr[i]);
  }

  sz = 0;
  er = fdf_readintarray(df,"g/e",&sz,err);
  quitOnError(*err,__LINE__,stderr);

  _DEBUGHERE_("g/e size %d",sz);
  for(i=0;i<sz;i++) {
    _DEBUGHERE_("g/e[%d] %ld",i,er[i]);
  }

  fdf_close(&df);

  _DEBUGHERE_("","");
}