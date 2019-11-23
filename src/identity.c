#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "memstack.h"
#include "mismatch.h"
#include "trimseq.h"
#include "alloc2d.h"
#include "bye.h"

#if 0
int vec_dim = 0;    /* maximum vector size */
int n_tau = 0;      /* number of tau values */
int n_theta = 0;    /* number of theta values */
double **poisson_tab=NULL; /* table of poisson distributions */
double *tau_tab=NULL;
double **eq_tab=NULL; /* table of equilibrium distributions */
double *theta_tab = NULL;
#endif

/** prototypes for internal use only */
int lookup(double value, double *vec, int vecdim);

#if 0
/****************************************************************
Call this before calling any of the other routines.  dim = size at
which vectors are to be allocated init_tau_tab is a vector of
length init_n_tau, containing all the values of tau that will
eventually be used.
****************************************************************/
void init_theory(int dim, double *init_tau_tab, int init_n_tau,
		 double *init_theta_tab, int init_n_theta)
{
  int i,j;

  vec_dim = dim;   /* all vectors will be allocated w/ this size */
  n_tau = init_n_tau;     /* # of tau values */
  n_theta = init_n_theta; /* # of theta values */

  tau_tab = (double *) mustalloc(n_tau * sizeof(double));
  theta_tab = (double *) mustalloc(n_theta * sizeof(double));
  poisson_tab = (double **) alloc2d((unsigned) n_tau, (unsigned) vec_dim,
				    sizeof(double));
  eq_tab = (double **) alloc2d((unsigned) n_theta, (unsigned) vec_dim,
			       sizeof(double));
  if(poisson_tab == NULL || eq_tab == NULL)
  {
    fflush(stdout);
    fprintf(stderr,"\ninit_theory: no memory\n");
    exit(1);
  }
  for(i=0; i<n_tau; i++)
    tau_tab[i] = init_tau_tab[i];
  qsort(tau_tab, (unsigned) n_tau, sizeof(double),
	(int (*)(const void*, const void*)) compar);

  for(i=0; i<n_tau; i++)
  {
    /* set poisson_tab[i] = exp(-tau)*tau^j/j! */
    poisson_tab[i][0] = exp(-tau_tab[i]);
    for(j=1; j<dim; j++)
      poisson_tab[i][j] = poisson_tab[i][j-1] * tau_tab[i] / j;
  }

  for(i=0; i<n_theta; i++)                            /* initialize vector */
    theta_tab[i] = init_theta_tab[i];
  qsort(theta_tab, (unsigned) n_theta, sizeof(double),
	(int (*)(const void*, const void*)) compar);  /* sort it */

  for(i=0; i<n_theta; i++)                        /* initialize eq dist's */
    geteq(eq_tab[i], vec_dim, theta_tab[i]);
}
/* Figure out how to initialize theory tables using POPHIST structure */
void hinit_theory(int msize, POPHIST *history)
{
  int i, dim;
  double *local_tau_tab, *local_theta_tab;
  POPHIST *ph;

  dim = 0;
  for(ph = history; ph != NULL; ph = ph->next)
    dim++;

  local_tau_tab = (double *) mustalloc((dim-1) * sizeof(double));
  local_theta_tab = (double *) mustalloc(dim * sizeof(double));

  i=0;
  for(ph = history; ph != NULL; ph = ph->next)
  {
    local_theta_tab[i] = ph->theta * ph->K;
    if(ph->next != NULL)
      local_tau_tab[i] = ph->tau;
    i++;
  }

  init_theory(msize, local_tau_tab, dim-1, local_theta_tab, dim);
  free(local_tau_tab);
  free(local_theta_tab);
}
#endif

/* un-initialize */
void free_theory(void)
{
  if(vec_dim == 0)
    return;
  vec_dim = n_tau = n_theta = 0;
  free2d((void **) poisson_tab);
  free2d((void **) eq_tab);
  free(tau_tab);
  poisson_tab = eq_tab = (double **) NULL;
  tau_tab = theta_tab = (double *) NULL;
}
/* lookup index of a value in a vector */
int lookup(double value, double *vec, int vecdim)
{
  int lo, mid, hi;

#ifndef NDEBUG
  assert(vecdim > 0);
  assert(vec != NULL);
#endif
#if 0
  printf("\nenter lookup %20.18g", value);
  for(mid=0; mid<vecdim && vec[mid] != value; mid++)
    ;
  if(mid==vecdim)
    printf("...can't find value");
  else
    printf("...=vec[%d]", mid);
  printf("\n%20.18g = value\ntable=", value);
  for(mid=0; mid<vecdim; mid++)
    printf("\n   %20.18g  (diff=%g)", vec[mid], value-vec[mid]);
  putchar('\n');
#endif
  lo = 0;
  hi = vecdim-1;
  
  while(lo<hi)
  {
    mid = lo + (hi-lo)/2;
    if(mid==lo)
    {
      if(value - vec[lo] > 0.5*(vec[hi] - vec[lo]))
	lo = hi;
      break;
    }
    if(value > vec[mid])
      lo = mid;
    else
      hi = mid;
  }
#if 1
  if(fabs(value - vec[lo]) > 1e-5*fabs(value))
  {
    printf("\nlookup error: |value-vec[lo]| > 1e-5*|value|");
    printf("\n%20.18g = vec[%d]", vec[lo], lo);
    printf("\n%20.18g = value\ntable=", value);
    for(mid=0; mid<vecdim; mid++)
      printf("\n   %20.18g  (diff=%g)", vec[mid], value-vec[mid]);
    putchar('\n');
    exit(1);
  }
#else
#ifndef NDEBUG  
  assert(value == vec[lo]);
#endif
#endif
  return(lo);
}


/****************************************************************
Equilibrium F distribution is geometric. theta= 2*mutation rate*pop size
****************************************************************/  
double *geteq(double *eq, int vecsize, double theta)
{
  double h;
  int i;

  eq[0] = 1.0/(theta + 1);
  
  h = 1.0 - eq[0];
  for(i=1; i<vecsize; i++)
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
double *nextf(double *f, double *f0, int vecsize, double theta, double tau)
{
  static double *diff;
  static int old_vecsize=0;
  double x, exp_factor;
  double *feq, *poi;
  int i, j;

#ifndef NDEBUG
  assert(tau >= 0.0);
  assert(vec_dim > 0);
  assert(n_tau > 0);
  assert(n_theta > 0);
#endif

  if(vecsize > old_vecsize)  /* allocate on 1st call */
  {
    if(vecsize > vec_dim)
    {
      fflush(stdout);
      fprintf(stderr,"\nError in nextf: vecsize=%d (max is %d)",
	      vecsize, vec_dim);
      exit(1);
    }
    if(old_vecsize > 0)
      free(diff);
    diff = (double *) mustalloc(vec_dim * sizeof(double));
    old_vecsize = vecsize;
  }
  feq = eq_tab[lookup(theta, theta_tab, n_theta)];

  if(tau==0.0)        /* resolve ambiguity for case when tau=theta=0 */
    exp_factor = 1.0;
  else
    exp_factor = exp(-tau/theta);

  /****************************************************************
  set diff[i] = f0[i] - feq[i]
  ****************************************************************/
  diff[0] = f0[0] - feq[0];
  for(i=1; i<vecsize; i++)
    diff[i] = f0[i] - feq[i];
  poi = poisson_tab[lookup(tau, tau_tab, n_tau)];

  for(i=0; i<vecsize; i++)
  {
    x = 0.0;
    for(j=0; j<=i; j++)
      x += poi[j] * diff[i-j];
    f[i] = feq[i] + x * exp_factor;
  }
  return(f);
}
/****************************************************************
Same as nextf(), except that the initial distribution is assumed to be
at equilibrium with theta = theta0.
****************************************************************/
double *getf(double *f, int vecsize, double theta0, double theta1, double tau)
{
#ifndef NDEBUG
  assert(n_theta > 0);
#endif
  return(nextf(f, eq_tab[lookup(theta0, theta_tab, n_theta)],
	       vecsize, theta1, tau));
}
/****************************************************************
Return mismatch distribution for a history list, ignoring population
structure. 
****************************************************************/
double *f_hist(POPHIST *ph, double *f, int dim)
{
#ifndef NDEBUG  
  assert(n_theta > 0);
#endif
  if(ph->next == NULL)
    return(eq_tab[lookup(ph->theta * ph->K, theta_tab, n_theta)]);
  /* else */
  return(nextf(f, f_hist(ph->next, f, dim),
	       dim, ph->theta * ph->K, ph->tau));
}
/****************************************************************
2-parameter mismatch distribution (theta1->infinity)
****************************************************************/
double *f_2param(double *f, int vecsize, double theta, double tau)
{
#if 1
  double *p, *eq;
#else
  double factor, *p;
#endif
  int i, j;

#ifndef NDEBUG
  assert(vecsize <= vec_dim);
  assert(n_tau > 0);
  assert(vec_dim > 0);
#endif

  p = poisson_tab[lookup(tau, tau_tab, n_tau)];
  eq = eq_tab[lookup(theta, theta_tab, n_theta)];
  /** calculate f **/
  for(i=0; i<vecsize; i++)
  {
    f[i] = 0.0;
#if 0
    factor = 1.0/(1.0+theta);
#endif
    for(j=i; j>=0; j--)
    {
#if 1      
      /* check to make sure that this does the same as the old code */
      f[i] += eq[i-j] * p[j];
#else
      f[i] += factor * p[j];
      factor *= theta/(1.0 + theta);
#endif
    }
  }
  return(f);
}
/****************************************************************
Add to final entry of mismatch distribution so that sum=1
****************************************************************/
double *fixsum(double *f, int vecsize)
{
  double cum=0.0;
  int i;

  /* make vector sum to unity by adding to terminal entry */
  for(i=0; i<vecsize; i++)
    cum += f[i];
  f[vecsize-1] += 1.0 - cum;
  return(f);
}

#if 0
main()
{
  int i,k=10, msize=50;
  double theta0, theta1;
  double tau;
  double *f;
  int t;

  f = (double *) mustalloc(msize * sizeof(double));

  for(;;)
  {
    printf("\nEnter theta0 theta1 tau: ");
    if(scanf("%f%f%f", &theta0, &theta1, &tau)==EOF)
      break;
    getf(f,k,theta0,theta1,tau);
    printf("\nf: ");
    for(i=0; i<k; i++)
      printf(" %f", f[i]);
  }
}
#endif
