#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "eprintf.h"
#include "mytypes.h"
#include "iscoales.h"
#include "mmtheory.h"
/*
 * The functions defined below calculate the theoretical mismatch
 * distribution.  For consistency, each function's arguments are in the
 * following order:
 *  foo(dest vec, src vec, vec dimension, parameters)
 */

/* get poisson distribution function */
double *getpoisson(double *poi, unsigned int dim, double tau)
{
  unsigned int j;

  if(dim==0 || poi==NULL)
    return(NULL);

  /* set poi[i] = exp(-tau)*tau^j/j! */
  poi[0] = exp(-tau);
  for(j=1; j<dim; j++)
    poi[j] = poi[j-1] * tau / j;

  return(poi);
}

/* get equilibrium mismatch distribution */
double *geteq(double *eq, unsigned int dim, double theta)
{
  double h;
  unsigned int i;

  if(dim==0 || eq==NULL)
    return(NULL);
  
  eq[0] = 1.0/(theta + 1);
  
  h = 1.0 - eq[0];
  for(i=1; i<dim; i++)
    eq[i] = h * eq[i-1];
  return(eq);
}

/****************************************************************
Calculate F distribution (f), given an arbitrary inititial
distribution (f0), and parameters theta (= 2*u*n, where u is mutation
rate and n the population size) and tau (= (2*u*t, where t is the
number of generations since the population was at its initial
distribution), and t (as just defined).  If t >= vecsize its value
will have no effect.  If f0 and f are the same vector, then the
original contents of f0 will be overwritten.

The formula used here is Eqn 4 of Rogers and Harpending (1992):

                                                    j
              ^                               i  tau            ^
  F (tau) =   F (tau) + exp(-tau(1+1/theta)) Sum --- (F   (0) - F   )
   i           i                             j=0  j!   i-j       i-j

or
              ^                          i                        ^
  F (tau) =   F (tau) + exp(-tau/theta) Sum Poi(j;tau) (F   (0) - F   )
   i           i                        j=0              i-j       i-j

where 
                   j
                tau
  Poi(j;tau) =  ---- exp(-tau)
                 j!

is the jth term of a Poisson distribution with parameter tau.
****************************************************************/
double *nextf(double *f, double *f0, unsigned int dim, double theta,
	      double tau)
{
  static double *diff;
  static double *feq, *poi;
  static unsigned int old_dim=0;
  double x, exp_factor;
  unsigned int i, j;
  char *progname = "nextf";

  if(dim==0 || f0==NULL || f==NULL)
    return(NULL);

  if(dim > old_dim) { /* allocate on 1st call */
    if(old_dim > 0) {
      free(diff);
      free(feq);
      free(poi);
    }
    diff = (double *) emalloc(progname, dim * sizeof(double));
    feq = (double *) emalloc(progname, dim * sizeof(double));
    poi = (double *) emalloc(progname, dim * sizeof(double));
    old_dim = dim;
  }
  geteq(feq, dim, theta);    /* get equilibrium distribution */
  getpoisson(poi, dim, tau); /* get poisson distribution */

  if(tau==0.0)        /* resolve ambiguity for case when tau=theta=0 */
    exp_factor = 1.0;
  else
    exp_factor = exp(-tau/theta);

  /*  set diff[i] = f0[i] - feq[i] */
  diff[0] = f0[0] - feq[0];
  for(i=1; i<dim; i++)
    diff[i] = f0[i] - feq[i];

  /* check this algorithm */
  for(i=0; i<dim; i++)  {
    x = 0.0;
    for(j=0; j<=i; j++)
      x += poi[j] * diff[i-j];
    f[i] = feq[i] + x * exp_factor;
  }
  return(f);
}
/* 
 * Same as nextf(), except that the initial distribution is assumed to be
 * at equilibrium with theta = theta0.
 */
double *getf(double *f, unsigned int dim, double theta0, double theta1,
	     double tau)
{
  static double *eq0=NULL;
  static unsigned int eq0_dim=0;
  char *progname = "getf";

  if(dim > eq0_dim) {  /* first time */
    if(eq0 == NULL)
      free(eq0);
    eq0 = (double *) emalloc(progname, dim * sizeof(double));
    eq0_dim = dim;
  }

  geteq(eq0, dim, theta0);

  return(nextf(f, eq0, dim, theta1, tau));
}
/*
 * Return mismatch distribution for a history list, ignoring population
 * structure.  This version allocates memory for temporary array and
 * then frees it.  It would be dangerous to use a static allocation
 * here because the recursive call to f_hist would have only one
 * vector rather than two.
 */
double *f_hist(double *f, unsigned int dim, POPHIST *ph)
{
  double *f0;
  char *progname = "f_hist";
  
  if(ph->next == NULL)
    return(geteq(f, dim, ph->theta * ph->K));
  /* else */
  f0 = (double *) emalloc(progname, dim * sizeof(double));
  f_hist(f0, dim, ph->next);
  nextf(f, f0, dim, ph->theta * ph->K, ph->tau);
  free(f0);
  return(f);
}

/* 2-parameter mismatch distribution (theta1->infinity) */
double *f_2param(double *f, unsigned int dim, double theta, double tau)
{
  static double *poi=NULL, *eq=NULL;
  static unsigned int old_dim=0;
  unsigned int i, j;
  char *progname = "f_2param";

  if(dim > old_dim) { /* 1st time */
    if(old_dim > 0) {
      free(poi);
      free(eq);
    }
    old_dim = dim;
    poi = (double *) emalloc(progname, old_dim * sizeof(double));
    eq = (double *) emalloc(progname, old_dim * sizeof(double));
  }

  geteq(eq, dim, theta);    /* get equilibrium distribution */
  getpoisson(poi, dim, tau); /* get poisson distribution */
  
  /** calculate f **/
  for(i=0; i<dim; i++)  {
    f[i] = 0.0;
    for(j=i; j>=0; j--)
      f[i] += eq[i-j] * poi[j];
  }
  return(f);
}

/* Add to final entry of vector so that sum=1 */
double *fixsum(double *f, unsigned int dim)
{
  double cum=0.0;
  unsigned int i;

  /* make vector sum to unity by adjusting terminal entry */
  for(i=0; i<dim; i++)
    cum += f[i];
  f[dim-1] += 1.0 - cum;
  return(f);
}

