#include "lollipop_lib.h"

int main()
{
  int i, nel, Lmin, Lmax, dl;
  double fsky, lnL;
  dataset *data;
  
  Lmin=2;
  Lmax=13;
  dl=2;
  fsky=0.77521502;

  char InFileName[256] = {0};
  sprintf( InFileName, "cross_100x143.dat");

  nel = (Lmax-Lmin+1)/dl;
  data = (dataset *) calloc( nel, sizeof(dataset));

  //Init datasets
  int status = Lollipop_Init( InFileName, Lmin, Lmax, dl, data);
  if( status !=0) {
    printf( "Error in Lollipop_Init\n");
    return(-1);
  }

  //read cl
  double tt, ee, bb, te;
  FILE *fp = fopen( "planck_beam_corbeam.bestfit_cl_zre8_r0.2.data", "r");
  unsigned int *ell  = (unsigned int *) malloc( (Lmax+1)*sizeof(unsigned int));
  double       *clee = (double       *) malloc( (Lmax+1)*sizeof(double      ));
  for( i=0; i<=Lmax; i++) {
    fscanf(fp,"%lf  %lf  %lf  %lf", &tt, &ee, &bb, &te);
    ell[i]  = i;
    clee[i] = ee*1e12;  //in muK2
  }
  fclose( fp);

  //compute Likelihood
  lnL = Lollipop_computeLikelihood( ell, NULL, NULL, clee, NULL, data, Lmin, Lmax, fsky, dl);
  printf( "Likelihood = %e\n", 2*lnL);

  return( 0);
}

