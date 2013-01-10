#include "fakedf.h"


void _read_meta(fdf *df,error **err) {
  char cuf[4096];
  char li[8192];
  char *li2,*li3;
  int ls;
  int nmax;
  FILE *ff;


  if (df->root[0]=='\0') {
    sprintf(cuf,"%s/%s",df->name,MDB);  
  } else {
    sprintf(cuf,"%s/%s/%s",df->root,df->name,MDB);  
  }
  
  ff = fopen_err(cuf,"r",err);
  forwardError(*err,__LINE__,);
  nmax = 10;

  df->metakey = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->metatype = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->metavalue = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->nmeta = 0;
    
  while(1) {
    int wh;
    ls = read_line(ff, li, 8192, err);
    testErrorRet(ls==-1,-1234,"bad line",*err,__LINE__,);  
    if (ls==1) {
      break;
    }
    
    if (df->nmeta==nmax) {
      resize_err(df->metakey, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      resize_err(df->metavalue, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      resize_err(df->metatype, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      nmax*=2;
    }

    for(wh=0;wh<ls;wh++) {
      if (li[wh]==' ') break;
    }
    li[wh]='\0';
    strcpy(df->metakey[df->nmeta],li);

    li2 = &(li[wh+1]);
    ls -= wh+1;
    for(wh=0;wh<ls;wh++) {
      if (li2[wh]==' ') break;
    }
    li2[wh]='\0';
    strcpy(df->metatype[df->nmeta],li2);

    li3 = &(li2[wh+1]);
    ls -= wh+1;
/*    for(wh=0;wh<ls;wh++) {
      if (li3[wh]==' ') break;
    }
    li3[wh]='\0';*/
    strcpy(df->metavalue[df->nmeta],li3);
    df->nmeta++;

  }
}


fdf * fdf_open_sub(char *path, char* sub,error **err) {
  fdf * df;
  struct stat buf;
  int erri;
  DIR *dd;
  struct dirent *dc;
  char *mpath;
  char mmpath[8192];
  int maxchild,maxdata;

  df = malloc_err(sizeof(fdf),err);
  
  forwardError(*err,__LINE__,NULL);
  df->root[0] = '\0';

  mpath = path;
  if (sub!=NULL && sub[0]!='\0') {
    strcpy(df->root,sub);  
    mpath = mmpath;
    sprintf(mpath,"%s/%s",sub,path);
  }

  erri = stat(mpath, &buf);
  testErrorRetVA(erri!=0,-1234,"cannot stat file (error %d)", *err,__LINE__,NULL,erri);
  testErrorRetVA(S_ISDIR(buf.st_mode)==0,-1234,"%s in not a directory",*err,__LINE__,NULL,path);

  strcpy(df->name,path);

  _read_meta(df,err);
  forwardError(*err,__LINE__,NULL);

  dd = opendir(mpath);
  testErrorRet(dd==NULL,-1234,"bad bad bad",*err,__LINE__,NULL);
  maxchild = 20;
  maxdata = 20;
  
  df->datakey = malloc_err(sizeof(kv)*maxdata,err);
  forwardError(*err,__LINE__,NULL);
  df->ndata=0;

  df->child = malloc_err(sizeof(void*)*maxchild,err);
  forwardError(*err,__LINE__,NULL);
  df->nchild=0;

  while(1) {
    dc = readdir(dd);
    if (dc==NULL) {
      break;
    }
    if (dc->d_type==DT_DIR && dc->d_name[0]!='.') {
      if (df->nchild==maxchild) {
        resize_err(df->child, sizeof(void*)*maxchild, sizeof(void*)*maxchild*2, 1, err);
        forwardError(*err,__LINE__,NULL);
        maxchild*=2;
      }
      df->child[df->nchild] = fdf_open_sub(dc->d_name,mpath,err);
      forwardError(*err,__LINE__,NULL);      
      df->nchild++;
    } else {
      if (dc->d_name[0]=='_') {
        continue;
      }
      if (df->ndata==maxdata) {
        resize_err(df->datakey, sizeof(void*)*maxdata, sizeof(void*)*maxdata*2, 1, err);
        forwardError(*err,__LINE__,NULL);
        maxdata*=2;
      }
      strcpy(df->datakey[df->ndata],dc->d_name);
      df->ndata++;
    }
  }
  return df;
}

fdf * fdf_open(char *path, error **err) {
  fdf *df;

  df = fdf_open_sub(path,NULL,err);
        forwardError(*err,__LINE__,NULL);
  return df;
}

fdf* fdf_tolast(fdf* df, char* key, char* kp, error **err) {
  int ls;
  int l0;
  int lk,i;
  fdf *cdf;

  l0 = 0;
  ls = 0;
  lk = strlen(key);
  cdf = df;

  while(1) {                
    if (key[lk]=='/') {     
      lk--;                 
    } else {                
      break;                
    }                       
  }             
  
  while(1) {
      
    if (key[ls]=='\0' || key[ls]=='/') {
      if (ls==l0) {
        l0++;
        ls=l0;
        continue;
      }
      memcpy(kp,&(key[l0]),sizeof(char)*(ls-l0));
      kp[ls-l0] = '\0';
      l0 = ls+1;
      ls = l0;
      
      // tester si c'est le dernier element
      if (ls>=lk) { //dernier element !
        return cdf;
      }
      // pas le dernier, deplace !
      for(i=0;i<cdf->nchild;i++) {
        if (strcmp(kp,((fdf*)cdf->child[i])->name)==0) {
            cdf = cdf->child[i];
            break;
        }
      }
    }
    ls++;
  }
}

int fdf_haskey(fdf *df, char *key, error **err) {
  char kp[256];
  int i;
  fdf * cdf;
  
  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      return 1;
    }
  }
  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      return 1;
    }
  }
  for(i=0;i<cdf->nchild;i++) {
    if (strcmp(kp,((fdf*) cdf->child[i])->name)==0) {
      return 1;
    }
  }
  return 0;
}

long fdf_readint(fdf *df, char *key, error **err) {
  char kp[256];
  int i;
  fdf * cdf;
  
  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"int")!=0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      return atol(cdf->metavalue[i]);          
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}


double fdf_readfloat(fdf *df, char *key, error **err) {
  char kp[256];
  int i;
  fdf * cdf;
  double res;

  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"str")==0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      testErrorRet(sscanf(cdf->metavalue[i],"%lg",&res)!=1,-1234,"gloups",*err,__LINE__,0);      
      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}

char* fdf_readstr(fdf *df, char *key, error **err) {
  char kp[256];
  int i;
  fdf * cdf;
  char* res;

  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);


  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"str")!=0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      res = malloc_err(sizeof(char)*(strlen(cdf->metavalue[i]+1)),err);
      forwardError(*err,__LINE__,NULL);
      strcpy(res,cdf->metavalue[i]);
      return res;
    }
  }
  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      FILE *f;
      struct stat st;
      char pth[8192];
      size_t size;
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      stat(pth, &st);
      size = st.st_size;
      res = malloc_err(size+1,err);
      forwardError(*err,__LINE__,NULL);
      f = fopen_err(pth,"r",err);
      fread(res, 1, size-1, f);
      res[size-1] = '\0';
      fclose(f);
      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}


fdf* fdf_openchild(fdf *df, char* key, error **err) {
  char kp[256];
  int i;
  fdf *cdf;

  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nchild;i++) {
    if (strcmp(kp,((fdf*)cdf->child[i])->name)==0) {
      return cdf->child[i];
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}

void fdf_close(fdf** pdf) {
  fdf *df,*cdf;
  int i;

  df = *pdf;
  if (df->root[0] == '\0') {
    free(df->metavalue);
    free(df->metakey);
    free(df->metatype);
    free(df->datakey);
    for(i=0;i<df->nchild;i++) {
      cdf = df->child[i];
      cdf->root[0] = '\0';
      fdf_close(&cdf);
    } 
  }
  *pdf = NULL;
}

void* fdf_readanyfits(char* path, int *sz, int typ, error **err) {
  void * res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero;
  int ploc;
  char ferrchar[80];
  int datatype;

  fitserr = 0;
  fits_open_data(&fitsptr, path, READONLY, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s) while opening %s",*err,__LINE__,NULL,fitserr,ferrchar,path);

  fitserr = 0;
  fits_get_img_type(fitsptr, &bitpix, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  testErrorRetVA(bitpix!=typ ,-1234,"bad type for file  %s",*err,__LINE__,NULL,path);

  fitserr = 0;
  fits_get_img_size(fitsptr, 1, &fsz, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  testErrorRetVA(fsz!=*sz && *sz>0 ,-1234,"bad size for file  %s (expected %d got %d)",*err,__LINE__,NULL,path,*sz,fsz);

  *sz = fsz;

  res = malloc_err((abs(typ)/8)*fsz,err);
  forwardError(*err,__LINE__,NULL);

  dzero = 0;
  fitserr = 0;
  datatype= TDOUBLE;
  if (typ>0) {
    datatype = TLONG;
  }
  fits_read_img(fitsptr, datatype, 1, fsz, &dzero, res, &ploc, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);

  fitserr = 0;
  fits_close_file(fitsptr, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while closing %s",*err,__LINE__,NULL,fitserr,ferrchar,path);

  return res;
 
}

double* fdf_readfloatarray(fdf *df, char *key, int* sz, error **err) {
  double *res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero; 
  int ploc;
  char kp[256];
  int i;
  fdf *cdf;

  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      char pth[8192];
      
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      res = fdf_readanyfits(pth,sz,-64,err);
      forwardError(*err,__LINE__,NULL);

      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}

long* fdf_readintarray(fdf *df, char *key, int* sz, error **err) {
  long *res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero;
  int ploc;
  char kp[256];
  int i;
  fdf *cdf;

  cdf = fdf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      char pth[8192];
      
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      res = fdf_readanyfits(pth,sz,64,err);
      forwardError(*err,__LINE__,NULL);
      
      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
}
