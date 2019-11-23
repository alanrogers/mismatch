#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "maxlike.h"

/* externals */
int first_time = 1;
TEST best;

/* record a test */
void record_test(double tau, double log10theta0, double growth, double pval)
{
  if(first_time || pval > best.pval)
  {
    best.tau = tau;
    best.log10theta0 = log10theta0;
    best.growth = growth;
    best.pval = pval;
  }
  first_time = 0;
}

/* get maximum likelihood estimates */
void getmaxlike(double *tau, double *theta0, double *theta1, double *pval)
{
  *tau = best.tau;
  *theta0 = pow(10.0, best.log10theta0);
  *theta1 = pow(10.0, best.growth) * *theta0;
  *pval = best.pval;
}

/* print maximum likelihood estimates in pictex format*/
void ptxmaxlike(FILE *fp)
{
  double tau, theta0, theta1, pval;

  getmaxlike(&tau, &theta0, &theta1, &pval);
  fprintf(fp, "\nMaximum likelihood estimates:");
  fprintf(fp, " $\\hat\\theta_0=%g$", theta0);
  fprintf(fp, ", $\\hat\\theta_1=%g$", theta1);
  fprintf(fp, ", $\\hat\\tau=%g$", tau);
  fprintf(fp, ", $\\hbox{$p$-val}=%g$", pval);
}

/* print maximum likelihood estimates in stdandard format*/
void prmaxlike(FILE *fp)
{
  double tau, theta0, theta1, pval;

  getmaxlike(&tau, &theta0, &theta1, &pval);
  fprintf(fp, "\nMaximum likelihood estimates:");
  fprintf(fp, " theta_0=%g", theta0);
  fprintf(fp, " theta_1=%g", theta1);
  fprintf(fp, " tau=%g", tau);
  fprintf(fp, " p-val=%g", tau);
}

