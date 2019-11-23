#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "memstack.h"
#include "mismatch.h"
#include "alloc2d.h"
#include "bye.h"
#include "chol.h"
#define MINITERATIONS 3
#undef TESTH0_VERBOSE
int testH0(SIMULATION *obs, double **x, int iterations, double eps,
	   FILE *fp, int msize, int dim_s, unsigned eflags, unsigned sflags)
{
  static double **c=NULL;  /* covariance matrix */
  static double *mean;     /* vector of means */
  static double *diff;
  static double xobs[MAXESTIMATES];
  double mean_seg, diff_seg;
  int    cholesky_err;
  double odist;
  double *v_mse;
  double *v_mae;
  double *v_roughness;
  double *v_seg;
  char *mlabel, *clabel;
  int i, j, count, nmahal, more_extreme_than_obs;

  if(iterations < MINITERATIONS)
    return(-3);


#ifdef TESTH0_VERBOSE
  printf("\n%c%21s", START_COMMENT, "testh0 entry");
  printf("\n%%testh0: obs->e.mse = %g", obs->e.mse);
#endif

  /* copy estimates from real data into xobs */
  i=j=0;
  if(sflags & THETA0) {
    if(eflags & THETA0)
      xobs[j++] = obs->e.theta0;
    else
      error("testH0: THETA0 is on in sflags but not in eflags");
  }
  if(sflags & THETA1) {
    if(eflags & THETA1)
      xobs[j++] = obs->e.theta1;
    else
      error("testH0: THETA1 is on in sflags but not in eflags");
  }
  if(sflags & TAU) {
    if(eflags & TAU)
      xobs[j++] = obs->e.tau;
    else
      error("testH0: TAU is on in sflags but not in eflags");
  }
  if(sflags & MSE) {
    if(eflags & MSE)
      xobs[j++] = obs->e.mse;
    else
      error("testH0: MSE is on in sflags but not in eflags");
  }
  if(sflags & MAE) {
    if(eflags & MAE)
      xobs[j++] = obs->e.mae;
    else
      error("testH0: MAE is on in sflags but not in eflags");
  }
  if(sflags & ROUGHNESS) {
    if(eflags & ROUGHNESS)
      xobs[j++] = obs->e.roughness;
    else
      error("testH0: ROUGHNESS is on in sflags but not in eflags");
  }
  if(sflags & SEG) {
    if(eflags & SEG)
      xobs[j++] = obs->e.seg;
    else
      error("testH0: SEG is on in sflags but not in eflags");
  }
  assert(j = countbits(sflags));

  /* count the variables to be used w/ Mahalanobis distance */
  nmahal = ((sflags&THETA0) !=0 ) + ((sflags&THETA1) != 0)
    + ((sflags&TAU) != 0);

  j = nmahal;
  if(sflags & MSE) 
    v_mse = x[j++];  /* vector of simulated mse values */

  if(sflags & MAE)   
    v_mae = x[j++];  /* vector of simulated mae values */

  if(sflags & ROUGHNESS)     
  v_roughness = x[j++];  /* vec of r'ghness vals */

  if(sflags & SEG) 
    v_seg = x[j++];  /* vec of segregating sites */

#ifdef TESTH0_VERBOSE
  printf("\n%cdim_s=%d nmahal=%d", START_COMMENT, dim_s, nmahal);
#endif  


#ifdef TESTH0_VERBOSE
  fprintf(fp, "\n%c %4s", START_COMMENT, "");
  fprintf(fp, "\n%c  %4s", START_COMMENT, "xobs");
  for(i=0; i<nmahal; i++)
    fprintf(fp, " %10.6f", xobs[i]);
#endif

  if(c==NULL && nmahal > 0) { /* allocate on first call */
    c = (double **) alloc2d(nmahal, nmahal, sizeof(double));
    if(c==NULL)
      error("testH0: can't allocate c");
    mean = (double *) mustalloc(nmahal * sizeof(double));
    diff = (double *) mustalloc(nmahal * sizeof(double));
  }

  /* get mean and covariance matrix */
  moments(x, nmahal, iterations, mean, c);  /* double mean & cov */
  mlabel = "Mean";
  clabel = "Covariance matrix";

  if(sflags & SEG) {
    /* get mean number of segregating sites */
    mean_seg = 0.0;
    for(i=0; i<iterations; i++)
      mean_seg += v_seg[i];
    mean_seg /= iterations;
    /* absolute diff between mean_seg and observed seg */
    diff_seg = fabs(obs->e.seg - mean_seg);
    if(fp!=NULL)
      fprintf(fp,"\n%c mean_seg=%f diff_seg=%f",
	      START_COMMENT, mean_seg, diff_seg);
  }

#if 0
  if(fp != NULL && nmahal > 0)  {
    fprintf(fp,"\n%cBefore rescaling c", START_COMMENT);
    fprintf(fp,"\n%c     ", START_COMMENT);
    for(i=0; i<nmahal; i++)
      fprintf(fp," %11s", lbl_s[i]);
    fprintf(fp,"\n%c%s:", START_COMMENT, mlabel);
    for(i=0; i<nmahal; i++)
      fprintf(fp," %11.8f", mean[i]);
    fprintf(fp,"\n%c%s:", START_COMMENT, clabel);
    prmat(c, nmahal, nmahal, fp);
  }
#endif

  /** rescale to improve accuracy **/
  for(i=0; i<nmahal; i++)  {
    for(j=0; j<nmahal; j++)
      c[i][j] /= mean[i]*mean[j];
  }

  /* observed difference vector */
  for(i=0; i<nmahal; i++)
    diff[i] = xobs[i]/mean[i] - 1.0;

#if 0
  /*******************************************************************
  This specifies the values of mean, c, and diff in order to make sure
  that the answer given by the version of chol and mahal used here is
  consistent with that given by other versions.  The answer should be
  odist = 84475
   *******************************************************************/
  assert(nmahal == 2);
  mean[0] = 3.73859890;
  mean[1] = 5.47871356;
  c[0][0] = 6.05230422 / (mean[0]*mean[0]);
  c[0][1] = c[1][0] = -0.19162061 / (mean[0]*mean[1]);
  c[1][1] = 1.68576091 / (mean[1]*mean[1]);
  diff[0] = -0.30481978;
  diff[1] = 68.76090210;
#endif
#if 0
  /*******************************************************************
  This specifies the values of mean, c, and diff in order to make sure
  that the answer given by the version of chol and mahal used here is
  consistent with that given by other versions.  The answer should be
  odist = 1.05178
   *******************************************************************/
  assert(nmahal == 2);
  mean[0] = 3.38834239;
  mean[1] = 5.59660784;
  c[0][0] = 0.52640522 ;
  c[0][1] = c[1][0] = -0.02725508 ;
  c[1][1] = 0.05051212  ;
  diff[0] = -0.23307042;
  diff[1] = 0.22788307;
#endif

#ifdef TESTH0_VERBOSE
  if(fp != NULL && nmahal > 0) {
    fprintf(fp,"\n%cAfter rescaling c", START_COMMENT);
    fprintf(fp,"\n%c     ", START_COMMENT);
    fprintf(fp,"\n%c%s:", START_COMMENT, "Obs");
    for(i=0; i<nmahal; i++)
      fprintf(fp," %11.8f", xobs[i]);
    fprintf(fp,"\n%c%s:", START_COMMENT, mlabel);
    for(i=0; i<nmahal; i++)
      fprintf(fp," %11.8f", mean[i]);
    fprintf(fp,"\n%c%s:", START_COMMENT, clabel);
    prmat(c, nmahal, nmahal, fp);

    fprintf(fp,"\n%cdiff: ", START_COMMENT);
    for(i=0; i<nmahal; i++)
      fprintf(fp, " %11.8f", diff[i]);
  }
#endif

  if(nmahal > 0) {
    cholesky_err = cholesky(c, nmahal); /* Cholesky factorization of c */
    if(cholesky_err)
      return(-1);   /* Cov mat not PD: no test is possible. */
    odist = mahal(diff, c, nmahal);  /* observed distance */
  }else
    odist = 1e6;
#ifdef TESTH0_VERBOSE
  fprintf(fp, "\n%codist = %g", START_COMMENT, odist);
#endif    

  count = 0;

  for(j=0; j<iterations; j++)  {
#ifdef TESTH0_VERBOSE
    fprintf(fp,"\n%cSimulation %d", START_COMMENT, j);
#endif
    for(i=0; i<nmahal; i++)
      diff[i] = x[i][j]/mean[i] - 1.0;

   /* Count simulations in tail: i.e. more extreme than obs */
    more_extreme_than_obs = 1;

    if(sflags & MSE) {
      /* 
       * Extreme cases fit better and have smaller mse.  If v_mse_2[j] >
       * obs->Mse_E, then the simulated value has a larger mse and is
       * less extreme.
       */
#ifdef TESTH0_VERBOSE
      fprintf(fp,"\n%c  MSE[%d]=%g obs=%g",
	      START_COMMENT, j, v_mse[j], obs->e.mse);
#endif
      if(more_extreme_than_obs && v_mse[j] > obs->e.mse)  {
#ifdef TESTH0_VERBOSE
	fprintf(fp,"\n%c  Failed MSE", START_COMMENT);
#endif
	more_extreme_than_obs = 0;
      }
    }

    if(sflags & MAE) {
      /*
       * Ditto for MAE.
       */
      if(more_extreme_than_obs && v_mae[j] > obs->e.mae) {
#ifdef TESTH0_VERBOSE
	fprintf(fp,"\n%c  Failed MAE", START_COMMENT);
#endif
	more_extreme_than_obs = 0;
      }
    }

    if(sflags & ROUGHNESS) {
      /* 
       * Extreme cases are very smooth and have low roughness values.
       * if v_roughness[j] > obs->Roughness_E, then the simulated value is
       * less smooth than the observed one and is therefore less
       * extreme.  
       */
      if(more_extreme_than_obs && v_roughness[j] > obs->e.roughness) {
#ifdef TESTH0_VERBOSE
	fprintf(fp,"\n%c  Failed Roughness", START_COMMENT);
#endif
	more_extreme_than_obs = 0;
      }
    }
    
    if(sflags & SEG) {
      /*
       * Extreme cases are farther from centroid than does the observed
       * number of segregating sites.  If fabs(v_seg[j]-mean_seg) < 
       * diff_seg, then case j is closer to the centroid and therefore
       * less extreme.
       */
      if(more_extreme_than_obs && fabs(v_seg[j]-mean_seg) < diff_seg)	{
#ifdef TESTH0_VERBOSE
	fprintf(fp,"\n%c  Failed SEG", START_COMMENT);
#endif
	more_extreme_than_obs = 0;
      }
    }
  
    /*
     * Extreme cases are farther from centroid, with correspondingly
     * larger Mahalanobis distance.  If mahal(diff, c, nmahal) < odist,
     * then the simulated case is closer to the centroid than the
     * observed case, so the simulated case is not more extreme.
     */
    if(more_extreme_than_obs && nmahal > 0 && mahal(diff, c, nmahal) < odist) {
#ifdef TESTH0_VERBOSE
      fprintf(fp,"\n%c  Failed Mahal", START_COMMENT);
#endif
      more_extreme_than_obs = 0;
    }
  
    if( more_extreme_than_obs )
      count++;
  }
  return( count );
}

void moments(double **x, int rows, int cols, 
		   double *mean, double **c)
{
  int i, j, k;
  double dx;

  /** mean **/
  for(i=0; i<rows; i++)
  {
    mean[i] = 0.0;
    for(k=0; k<cols; k++)
      mean[i] += x[i][k];
    mean[i] /= cols;
  }

  for(i=0; i<rows; i++)
  {
    /** variances **/
    c[i][i] = 0.0;     
    for(k=0; k<cols; k++)
    {
      dx = x[i][k] - mean[i];
      c[i][i] += dx*dx;
    }
    c[i][i] /= cols - 1.0;

    /** covariances **/
    for(j=0; j<i; j++)
    {
      c[i][j] = 0.0;
      for(k=0; k<cols; k++)
	c[i][j] += (x[i][k] - mean[i])*(x[j][k] - mean[j]);
      c[i][j] /= cols - 1.0;
      c[j][i] = c[i][j];
    }
  }
}

