// this fakes hdf5 API with a bunch of directories along with fits and text files

#include "errorlist.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "fitsio.h"

#define MDB "_mdb"

typedef char kv[256];

typedef struct {
  char name[2048];
  char root[2048];
  int nmeta,ndata,nchild;
  kv* metakey;
  kv* metavalue;
  kv* metatype;
  kv* datakey;
  void **child;
} fdf;

fdf * fdf_open(char *path, error **err);
fdf* fdf_openchild(fdf *df, char* ipath, error **err);
int fdf_haskey(fdf *df, char *key, error **err);
long fdf_readint(fdf *df, char *key, error **err);
double fdf_readfloat(fdf *df, char *key, error **err);
char* fdf_readstr(fdf *df, char *key, error **err);
long* fdf_readintarray(fdf *df, char *key, int* sz, error **err);
double* fdf_readfloatarray(fdf *df, char *key, int* sz, error **err);

void fdf_close(fdf** pdf);

