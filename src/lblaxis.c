#include <stdio.h>
#include <math.h>
#include "mytypes.h"

/**** prototypes  ****/
void lblaxis(double *min, double *max, int *ntic, char **ticmark,
	     double *ticvalue, int logscale);
void printwidth(double x, int *len, int *frac);
void pictex_ticks(FILE *fp, int nticks, char **tickmark, double *tickvalue);
double tryticscale(double target_width, double rval, int exp10);

#define MAXXCESS 0.5
#define NRVAL 4
double roundval[NRVAL] = {1.0, 2.0, 2.5, 5.0};

/****************************************************************
Do calculations for pretty axis labels.

On entry: askmin, askmax, and askntic give the requested minimum tick
mark, maximum tick mark, and number of tick marks.  Note that these
are all pointers since they will be modified by the function.

On return:

 *askmin contains the actual minimum tick mark
 *askmax contains the actual maximum tick mark
 *askntic contains the number of tick marks
 tic  is a vector containing the tick marks
 fmt  is a character string containing the printf format that should be
      used to print them.  In other words, print the i'th tick mark
      with:  "printf(fmt, tic[i])".
****************************************************************/
void lblaxis(double *askmin, double *askmax, int *askntic, char **ticmark,
	     double *ticvalue, int logscale)
{
  double range, w, x;
  double bestrval;
  int try_exp10, exp10, bestexp10;
  double bestfit, fit;
  double max, min, oldmax, oldmin;
  int i, n, len, frac;
  char fmt[50];

  max = *askmax;
  min = *askmin;
  n   = *askntic;

  do{
    oldmax = max;
    oldmin = min;
    range = max - min;
    w = range / (n - 1.0);
    /****************************************************************
    Add 0.1 on next line to make sure roundoff error doesn't move us
    to next smaller int.
    ****************************************************************/
    bestexp10 = exp10 = floor(log10(w)) + 0.1;  
    bestrval = roundval[0];
    bestfit = tryticscale(w, roundval[0], exp10-1);
    for(try_exp10 = exp10-1; try_exp10 <= exp10+1; try_exp10++)
    {
      for(i=0; i<NRVAL; i++)
      {
	fit = tryticscale(w, roundval[i], try_exp10);
	if(fit < bestfit)
	{
	  bestfit = fit;
	  bestrval = roundval[i];
	  bestexp10 = try_exp10;
	}
      }
    }
    /* width of interval is roundval[i0] * 10^exp10 */
    w = bestrval * pow(10.0, (double) bestexp10);
    
    /* calculate new min */
    min = ceil(*askmin/w) * w;
#if 0    
    if(min - *askmin > MAXXCESS*(*askmax - *askmin))
      min = (ceil(*askmin/w)-1.0) * w;     /* new min */
#endif    
    
    /* calculate new max */
    max = floor(*askmax/w) * w;
#if 0    
    if(*askmax - max > MAXXCESS*(*askmax - *askmin))
      max = (floor(*askmax/w)+1.0) * w;
    }while(oldmax != max || oldmin != min);
#endif    
  }while(0);
    
  *askntic = n = (max - min)/w + 1.5; /* rounds to nearest int */

  /* fill tick-value vector */
  for(i=0; i<n; i++)
    ticvalue[i] = min + i*w;

  *askmin = ticvalue[0];
  *askmax = ticvalue[n-1];

  /* fill tick-mark vector */
  for(i=0; i<n; i++)
  {
    if(logscale)
      x = pow(10.0,ticvalue[i]);
    else
      x = ticvalue[i];
    printwidth(x, &len, &frac);
    sprintf(fmt,"%%%d.%df", len, frac);
    sprintf(ticmark[i], fmt, x);
  }
}

double tryticscale(double target_width, double rval, int exp10)
{
  return(fabs( target_width - rval * pow(10.0, (double) exp10)));
}
void printwidth(double x, int *len, int *frac)
{
  int left, right, sign, decimal=0;
  double intpart, fracpart;

  if(x < 0.0)
  {
    sign = -1;
    x = -x;
  }else
    sign = 1;

  intpart = floor(x);
  fracpart = x - intpart;
  if(intpart == 0.0 || intpart==1.0)
    left = 1;                         /* space for leading digit */
  else
    left = ceil(log10(intpart));      /* space for integer part */
  if(fracpart == 0.0)
    right = 0;
  else
    right = ceil(log10(1.0/fracpart));  /* space for fractional part */
  if(sign == -1)
    left++;                           /* space for minus sign */
  if(right > 0)
    decimal = 1;                      /* space for decimal point */
    
  *len = left + decimal + right;      /* total printing width */
  *frac = right;                      /* width of fractional part */
}
/** Pictex output **/
void pictex_ticks(FILE *fp, int nticks, char **tickmark, double *tickvalue)
{
  int i;
  
  fprintf(fp,"\n   ticks withvalues");
  for(i=0; i<nticks; i++)
    fprintf(fp, " %s", tickmark[i]);
  fprintf(fp," /\n     at");
  for(i=0; i<nticks; i++)
    fprintf(fp, " %f", tickvalue[i]);
  fprintf(fp," / /");
}
#if 0
void main(void)
{
  double min, max, tic[100];
  int i, ntic;
  char fmt[50];
  for(;;)
  {
    fprintf(stderr,"\nEnter min max ntic:");
    fscanf(stdin,"%f%f%d", &min, &max, &ntic);
    fprintf(stderr,"\nEcho input: min=%f max=%f ntic=%d", min, max, ntic);
    lblaxis(&min, &max, &ntic, tic, fmt);
    for(i=0; i<ntic; i++)
      fprintf(stderr,"\n %g", tic[i]);
    fprintf(stderr,"\nFmt=%s", fmt);
  }
}
#endif
