#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct {
  double l, cl, na, nb;
} dataset;

int Lollipop_Init( char *InFileName, int Lmin, int Lmax, int dl, dataset *data);
double Lollipop_H2( double y);
double Lollipop_H3( double y);
double Lollipop_EdgeWorthSeries( int dl, double fsky, double cl0, dataset data);
double Lollipop_computeLikelihood(const unsigned int *l,
          double *cltt,
          double *clte,
          double *clee,
          double *clbb,
          dataset *data, 
          int Lmin, int Lmax, double fsky, int dl);



