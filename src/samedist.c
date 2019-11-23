#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "alloc2d.h"
#include "unirand.h"
#include "samedist.h"
#include "bye.h"
/* DMAT simplifies use of lower triangular matrix dmat */
#define DMAT(i,j) ( (i)>(j) ? dmat[i][j] : dmat[j][i] )

/****************************************************************
Test the hypothesis that two samples are drawn from the same
distribution.  

On input: 
  s1   is an n2 X dim matrix of doubles.  The i'th row is the
       i'th observation from distribution 1.
  s2   is an n2 X dim matrix, whose i'th row contains the i'th
       observation from distribution 2.
  n1   is the number of observations from distribution 1.
  n2   is the number of observations from distribution 2.
  dim  is the dimension of each observation.
  nrand
       number of randomizations to be performed
  dist is a function whose three arguments are
         obs1 a vector of dim doubles
         obs2 another such vector
	 dim  the length of each vector
       The function dist should return the distance between obs1
       and obs2.
       
On return from samedist:
  The function returns the p-value of the null hypothesis, that is, the
hypothesis that the two samples are drawn from the same distribution.
If p < 0.05 then the null hypothesis can be rejected at the 0.05
significance level.

Algorithm: Calculate distance between each pair of observations in the
combined data set, and rescale these as proximities, where the
proximity of i and j is 0 if the distance between them exceeds 'chop',
and equals exp(-scale*distance) if distance<=chop.  I set chop equal
to the mean w/i sample distance and scale = 4/chop.  Next calculate

 diff_mean = (mean w/i sample proximity) - (mean btw sample proximity)

Initialize the variable 'tail' and then repeat the following  
nrand times:

  1. Form a vector ndx that is a random permutation of 0..(n1+n2-1).

  2. Use ndx to shuffle the observations and then calculate
  r_diff_mean, which is defined like diff_mean but with the randomly
  shuffled observations.

Return the fraction of randomizations in which r_diff_mean >= diff_mean.

Under the null hypothesis, the assignment of observations to groups is
random, so r_diff_mean and diff_mean are drawn from the same
distribution.

Alan Rogers  alan.rogers@anthro.utah.edu     7 February 1997
****************************************************************/
double samedist(double **s1, int n1, double **s2, int n2, int dim,
	       int nrand,
	       double (*d)(double *obs1, double *obs2, int dim) )
{
  static int n=0;
  static int *neighbor, *ndx;
  static double **dmat;
  double chop, meanwi, meanbtw, diff_mean, r_diff_mean;
  double scale;
  int i,j, tail, iteration;

#if 0
  printf("\nsamedist: sample1 (%d):", n1);
  for(i=0; i<n1; i++)
  {
    putchar('\n');
    for(j=0;j<dim;j++)
      printf(" %5g", s1[i][j]);
  }

  printf("\nsamedist: sample2 (%d):", n2);
  for(i=0; i<n2; i++)
  {
    putchar('\n');
    for(j=0;j<dim;j++)
      printf(" %5g", s2[i][j]);
  }
#endif

  if(n < n1+n2)  /* allocate arrays */
  {
    if(n>0)
    {
      free(neighbor);
      free2d((void **) dmat);
    }
    n = n1+n2;
    neighbor = (int *) malloc(n * sizeof(int));
    ndx = (int *) malloc(n * sizeof(int));
    dmat = (double **) alloclt(n, sizeof(double));
    if(neighbor==NULL || ndx==NULL || dmat==NULL )
      error("memory");
  }
  n = n1+n2;

  /* fill distance matrix */
  chop = 0.0;
  for(i=1; i<n1; i++)
    for(j=0; j<i; j++)
      chop += dmat[i][j] = (*d)(s1[i], s1[j], dim);

  for(i=0; i<n2; i++)
    for(j=0; j<n1; j++)
      dmat[n1+i][j]  = (*d)(s2[i], s1[j], dim);

  for(i=0; i<n2; i++)
    for(j=0; j<i; j++)
      chop += dmat[n1+i][n1+j] = (*d)(s2[i], s2[j], dim);

  chop /= n1*(n1-1)/2 + n2*(n2-1)/2;  /* chop = mean dist wi */
  /****************************************************************
  Set scale so that exp(-scale*x) is small when x > chop
  ****************************************************************/
  scale = 4.0/chop;

  /* convert distance matrix into proximity matrix */
  for(i=1; i<n; i++)
    for(j=0; j<i; j++)
    {
      if(dmat[i][j] >= chop)
	dmat[i][j] = 0.0;
      else
	dmat[i][j] = exp(-scale*dmat[i][j]);
    }

  /* calculate w/i and btw sample sums */
  meanwi = meanbtw = 0.0;
  for(i=1; i<n1; i++)
    for(j=0; j<i; j++)
      meanwi += dmat[i][j];

  for(i=0; i<n2; i++)
    for(j=0; j<n1; j++)
      meanbtw += dmat[n1+i][j];  

  for(i=0; i<n2; i++)
    for(j=0; j<i; j++)
      meanwi += dmat[n1+i][n1+j];

  meanwi /= n1*(n1-1)/2 + n2*(n2-1)/2;
  meanbtw /= n1*n2;
  diff_mean = meanwi - meanbtw;

  for(tail = iteration = 0; iteration < nrand; iteration++)
  {
    (void) randperm(ndx, n);  /* random permutation of 0..(n-1) */

    meanwi = meanbtw = 0.0;
    for(i=1; i<n1; i++)
      for(j=0; j<i; j++)
	meanwi += dmat[ndx[i]][ndx[j]];

    for(i=0; i<n2; i++)
      for(j=0; j<n1; j++)
	meanbtw += dmat[ndx[n1+i]][ndx[j]];

    for(i=0; i<n2; i++)
      for(j=0; j<i; j++)
	meanwi += dmat[ndx[n1+i]][ndx[n1+j]];

    meanwi /= n1*(n1-1)/2 + n2*(n2-1)/2;
    meanbtw /= n1*n2;
    r_diff_mean = meanwi - meanbtw;
    if(r_diff_mean >= diff_mean)
      tail += 1;
  }
  return((double) tail/ (double) nrand);
}

double dist(double *obs1, double *obs2, int dim)
{
  double sum=0.0;
  int i;

  for(i=0; i<dim; i++)
  {
    sum += fabs( obs1[i] - obs2[i] );
  }
  return(sum);
}

#if 0
/**************************************************************** 
Driver routine for testing.  For real applications, turn this off
by changing '#if 1' to '#if 0' on the line above.
****************************************************************/
void main(void)
{
  double **s1, **s2, pvalue;
  int i, n1=200, n2=200, dim=1;
  int nreject, rep, nrep = 10;
  double mean, sd;

  fputs("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", stdout);

  initrand(0);

  s1 = (double **) alloc2d(n1, dim, sizeof(double));
  s2 = (double **) alloc2d(n2, dim, sizeof(double));
  if(s1 == NULL || s2 == NULL)
    error("memory");

  mean = 0.75;
  sd = 1.0;

  printf("\nsample 1: Normal(0,1) n=%d", n1);
  printf("\nsample 2: 0.5*Normal(%g,%g) + 0.5*Normal(-%g,%g) n=%d",
	 mean, sd, mean, sd, n2);
  printf("\ndim=%d", dim);

  for(nreject=rep=0; rep<nrep; rep++)
  {
    for(i=0; i<n1; i++)
      s1[i][0] = rnorm();

    /* sample from normal(mean, sd) */
    for(i=0; i<n2; i++)
    {
      if(uni() > 0.5)
	s2[i][0] = mean + sd*rnorm();
      else
	s2[i][0] = -mean + sd*rnorm();
    }
    pvalue = samedist(s1, n1, s2, n2, dim, 500, dist);

    printf("\n%3d: pvalue=%f", rep+1, pvalue);
    if(pvalue <= 0.05)
      nreject++;
    printf(" rejected: %d/%d", nreject, rep+1);
  }
  putchar('\n');
}
#endif

