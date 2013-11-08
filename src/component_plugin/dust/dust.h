#include "clik_parametric.h"
#include "clik_parametric_addon.h"
#define PRM_NU0 143.

void gal_beta_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void gal_T_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void gal_alpha_non_thermal_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
