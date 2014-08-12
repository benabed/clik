#include "clik.h"
#include "clik_helper.h"
#include "lollipop_lib.h"

typedef struct {
  dataset* data;
  int dl,nel,lmin,lmax;
  double fsky;
  unsigned int *ell;
  double *cltt;
  double *cls[6];
  int has_cl[6];
} lollipop_payload;

void free_clik_lollipop(void **plol) {
  lollipop_payload *lol;

  lol = *plol;
  free(lol->ell);
  free(lol->data);
  free(lol->cls[0]);
  free(lol);
  *plol = NULL;
}

double clik_lollipop_lkl(void* vlol, double* pars, error **err) {
  lollipop_payload *lol;
  int i,cli,tot;
  double res;

  lol = vlol;
  tot = 0;
  for(cli=0;cli<6;cli++) {
    if (lol->has_cl[cli]==0) {
      continue;
    }
    for(i=lol->lmin;i<lol->lmax+1;i++) {
      lol->cls[cli][i] = pars[tot];
      tot++;
    }
  }
  res = Lollipop_computeLikelihood( lol->ell, lol->cls[0], lol->cls[3], lol->cls[1], lol->cls[2], lol->data, lol->lmin, lol->lmax, lol->fsky, lol->dl);
  return -res;
}

cmblkl* clik_lollipop_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cmblkl *cing;
  lollipop_payload *lol;
  int i,status;
  char directory_name[4096],pwd[4096];
  double *cltt;

  lol = malloc_err(sizeof(lollipop_payload),err);
  forwardError(*err,__LINE__,NULL);

  lol->dl = cldf_readint(df,"dl",err);
  forwardError(*err,__LINE__,NULL);
  lol->nel = cldf_readint(df,"nel",err);
  forwardError(*err,__LINE__,NULL);
  lol->fsky = cldf_readfloat(df,"fsky",err);
  forwardError(*err,__LINE__,NULL);
  lol->lmin = ell[0];
  lol->lmax = ell[nell-1];

  for(i=0;i<6;i++) {
    lol->has_cl[i] = has_cl[i];
  }

  lol->ell = malloc_err(sizeof(unsigned int)*(lol->lmax+1),err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<lol->lmax+1;i++) {
    lol->ell[i] = i;
  }
  
  cltt = malloc_err(sizeof(double)*(lol->lmax+1)*6,err);
  forwardError(*err,__LINE__,NULL);

  memset(cltt,0,sizeof(double)*(lol->lmax+1)*6);
  for(i=0;i<6;i++) {
    lol->cls[i] = cltt + (lol->lmax+1)*i;
  }

  lol->data = calloc_err( lol->nel, sizeof(dataset),err);

  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  status = Lollipop_Init( "lollipop_data.dat", lol->lmin, lol->lmax, lol->dl, lol->data);
  testErrorRetVA(status!=0,-324324,"bad status after init (reported %d)",*err,__LINE__,NULL,status);

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  cing = init_cmblkl(lol, &clik_lollipop_lkl, 
                     &free_clik_lollipop,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  
  return cing;
}

