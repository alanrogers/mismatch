#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "memstack.h"
#include "mismatch.h"
#include "estimate.h"
#include "eprintf.h"
#include "mmtheory.h"
/****************************************************************
If theory_f==NULL then Calculate MSE using estimated parameters.
Otherwise, use distribution in theory_f.
****************************************************************/
void estimate(SIMULATION *s, int msize, double *theory_f, int *match,
	      unsigned eflags)
{
  double theta[2], *ff;
  int small_dim;
  static int old_msize = 0;
  static double *f;
  char *progname = "estimate";
  
#undef ESTIMATE_VERBOSE

#ifdef ESTIMATE_VERBOSE
  printf("\n%cestimate entry: msize=%d theory_f=%u old_msize=%d f=%u",
	 START_COMMENT, msize, (unsigned) theory_f, old_msize, (unsigned) f);
#endif

  if(theory_f == NULL && msize != old_msize)
  {
    if(old_msize > 0)
      free(f);
    old_msize = msize;
    f = (double *) emalloc(progname, msize * sizeof(double) );
    ff = f;
#ifdef ESTIMATE_VERBOSE
    printf("\n%cestimate: reallocating f, old_msize set = %d",
	   START_COMMENT, old_msize);
#endif
  }else
    ff = theory_f;
  
  getcumulants(match, msize, s->e.cumulant);

  if( eflags & (THETA0 | THETA1 | TAU)
      || ( eflags & (MSE|MAE) && theory_f==NULL))
  {
#ifdef ESTIMATE_VERBOSE
    printf("\n%c estimate: calling mom2", START_COMMENT);
#endif
    /* estimate theta0, theta1, tau */
    mom2(&(s->e.theta0), &(s->e.theta1), &(s->e.tau),
	 s->e.cumulant, match, msize);
    if( ( eflags & (MSE|MAE) ) && theory_f==NULL  )
    {
#ifdef ESTIMATE_VERBOSE
      printf("\n%c estimate: calculating mse from parameter estimates",
	     START_COMMENT);
#endif
      small_dim = veclength(match, msize);  /* trim zeroes from end */
      theta[0] = s->e.theta0;
      theta[1] = s->e.theta1;
      getf(f, msize, s->e.theta0, s->e.theta1, s->e.tau);
      if( eflags & MSE )
	s->e.mse = get_mse(f, match, msize);  /* set MSE */
      if( eflags & MAE )
	s->e.mae = get_mae(f, match, small_dim);  /* set MAE */
    }
  }
  if( (eflags & MSE) && theory_f!=NULL)   /* set MSE using theory_f */
  {
#if 0
    printf("\n%c estimate: using theory_f for MSE\n%%  theory_f:",
	   START_COMMENT);
    for(i=0; i<msize; i++)
      printf(" %g", theory_f[i]);
#endif
    s->e.mse = get_mse(theory_f, match, msize);
  }
  if( (eflags & MAE) && theory_f!=NULL)   /* set MAE using theory_f */
  {
#ifdef ESTIMATE_VERBOSE
    printf("\n%% estimate: using theory_f for MAE\n%%  f:");
    for(i=0; i<msize; i++)
      printf(" %g", f[i]);
#endif
    s->e.mae = get_mae(theory_f, match, msize);
  }
  if( eflags & ROUGHNESS )
    s->e.roughness = getroughness(match, msize);
}
/****************************************************************
  2-parameter method of moments estimator. 
****************************************************************/   
void mom2(double *theta0, double *theta1, double *tau,
	  double *k, int *h, int msize)
{
  double sum;
  int i;

#ifndef NDEBUG
  assert(k[0] >= 0.0);
  assert(k[1] >= 0.0);
#endif

  if(k[0] > k[1])  /* mean > variance */
    *theta0 = 0.0;
  else
    *theta0 = sqrt(k[1] - k[0]);

  for(sum=0.0,i=0; i<msize; i++)
    sum += h[i];
  
  *tau = k[0] - *theta0;
  if(*tau < 0.0)
    *tau = 0.0;
  *theta1 = (sum - h[0])/h[0];
}

/* Harpending's roughness statistic, r, is the sum of squared differences */
double getroughness(int *h, int msize)
{
  int i;
  double sum = 0.0, sumsq=0.0, dif;

  if(h == NULL || msize < 2)
    return(0.0);           /* roughness undefined */
  for(i=0; i<msize; i++)
    sum += h[i];           /* sum = sum of mismatch distribution */
  for(i=0; i<msize-1; i++)
  {
    dif = (h[i+1] - h[i])/sum;  /* difference between adjacent entries */
    sumsq += dif*dif;           /* sum of squares */
  }
  return(sumsq);         /* return sum of squared diffs */
}

/* mean squared error */
double get_mse(double *f, int *h, int msize)
{
  double rsumh, diff, sumdiff=0.0;
  int i, sumh=0;

  for(i=0; i<msize; i++)
    sumh += h[i];

#ifndef NDEBUG
  assert(sumh != 0);
#endif
  rsumh = (double) sumh;    /* convert int to double */

  for(i=0; i<msize; i++)
  {
    diff = f[i] - h[i]/rsumh;
    sumdiff += diff*diff;
  }
#if 0  
  if(sumdiff == 0)
  {
    printf("\nget_mse found 0 diffs.  sumh=%g:", sumh);
    for(i=0; i<msize; i++)
      printf("\n  %d: f=%f h=%f", i, f[i], h[i]/sumh);
  }
#endif  
  return(sumdiff / msize);
}

/****************************************************************
mean absolute error
f = theoretical mismatch distribution
h = empirical mismatch distribution (integers)
dim = length of vectors
****************************************************************/
double get_mae(double *f, int *h, int dim)
{
  double rsumh, cumf=0.0, sumdiff=0.0;
  int i, sumh=0, hsize;
#if 0
  extern double obs_mae;
#endif

  for(hsize = dim; h[hsize-1]==0; hsize--)  /* exclude terminal zeroes */
    ;
  for(i=0; i<hsize; i++)
  {
    sumh += h[i];
    cumf += f[i];
  }

#ifndef NDEBUG
  assert(sumh != 0);
#endif
  rsumh = (double) sumh;    /* convert int to double */

  for(i=0; i<hsize; i++)
    sumdiff += fabs(rsumh*f[i] - h[i]);
  sumdiff /= rsumh;

  if(hsize < dim)
    sumdiff += 1.0 - cumf;
#ifndef NDEBUG
  if(hsize==dim && fabs(1.0-cumf) > 0.001)
  {
    printf("\nNumerical error in get_mae: hsize=dim=%d but cumf=%g (!=1)",
	   dim, cumf);
    printf("\nf=");
    for(i=0; i<dim; i++)
      printf(" %g", f[i]);
    putchar('\n');
    exit(1);
  }
#endif
#if 0
  if(sumdiff / dim >= obs_mae)
  {
    printf("\nMAE = %g >= %g = obs_MAE", sumdiff/dim, obs_mae);
    printf("\nf=");
    for(i=0; i<hsize; i++)
      printf(" %6.4g",f[i]);
    printf("\nh=");
    for(i=0; i<hsize; i++)
      printf(" %6.4g", h[i]/rsumh);
  }
#endif
  return(sumdiff / dim);
}
