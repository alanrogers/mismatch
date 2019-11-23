#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "charplot.h"
#include "bye.h"
#define OW 78			/* output width */
#define OH 24			/* output height */
#define MAXTIC 25		/* maximum # of tics */


/*--------------------CHARPLOT-------------------------------------
Plots on terminals w/o graphics capabilities.
On entry: 
     fp is pointer to output file
      n is length of data vectors
      x is a double vector of x-axis values
      y is a double vector of y-axis values
Function writes plot on file fp.
------------------------------------------------------------------------*/

void charplot(FILE *fp, int n, double *x, double *y)
{
  int     xmin, xmax, ymin, ymax, i, j, yzero, xzero,
          xm, ym;
  char    map[OH + 1][OW + 1];

  double   xtic[MAXTIC];
  int     nxtic = 0, nytic = 0;

  /* Initialize parameters -------------------------- */
  xm = ym = 0;

  xmin = 1 + xm;
  xmax = OW - 2 - xm;
  ymin = 1 + ym;
  ymax = OH - 1 - ym;

  for(i=0; i<MAXTIC; i++)
    xtic[i] = 10*i;
  nxtic = MAXTIC;
  nytic = 0;

  /* Adjust x ,y and tics to fit on map ------------- */
  xzero = adjust(x, n, xmin, xmax, xtic, &nxtic);
  yzero = adjust(y, n, ymin, ymax, NULL, &nytic);

  /* create plot matrix--------------------------------- */
  for (i = 0; i <= OH; i++)	/* initialize map */
    for (j = 0; j <= OW; j++)
      map[i][j] = ' ';
  for (i = 1; i < OH; i++)	/* draw side lines */
    map[i][0] = map[i][OW - 1] = '|';
  for (i = 0; i < OW; i++)	/* draw top & bottom lines */
    map[0][i] = map[OH][i] = '-';
#if 0
  for (i = 0; i < nytic; i++)	/* tics on y axis */
  {
    if (ytic[i] >= 0 && ytic[i] < OH)
      map[(int) ytic[i]][0] = map[(int) ytic[i]][OW - 1] = '+';
  }
#endif
  for (i = 0; i < nxtic; i++)	/* tics on x axis */
  {
    if (xtic[i] >= 0 && xtic[i] <= OW - 1)
      map[0][(int) xtic[i]] = map[OH][(int) xtic[i]] = '+';
  }
  for (i = 0; i < n; i++)	/* data points */
    map[(int) (y[i]+0.5)][(int) (x[i]+0.5)] = '*';
  for (i = 0; i <= OH; i++)	/* terminate lines w/ null */
    map[i][OW] = '\0';

  /*---------Print map----------------------------------------*/
  for (i = OH; i >= 0; i--)
  {
    putc('\n', fp);
    fputs(map[i], fp);
  }
}

/*
 * Rescale values of x so that they range from min to max.  Use same
 * rescaling on values of xtic.
 */
int adjust(double *x, int n, int min, int max, double *xtic, int *nxtic)
{
  double   xmin, xmax, scalefactor;
  int     i, zero;

  xmin = xmax = x[0];
  for (i = 1; i < n; i++)	/* find max and min of x */
  {
    if (x[i] > xmax)
      xmax = x[i];
    if (x[i] < xmin)
      xmin = x[i];
  }
  if (xmax == xmin)
    scalefactor = 1;
  else
    scalefactor = (max - min) / (xmax - xmin);

  for (i = 0; i < n; i++)
    x[i] = round((x[i] - xmin) * scalefactor) + min;
  if(xtic != NULL)
  {
    for (i = 0; i < *nxtic; i++)
    {
      if(xtic[i] > xmax)
      {
	*nxtic = i;
	break;
      }
      xtic[i] = round((xtic[i] - xmin) * scalefactor) + min;
    }
  }
  if (xmin > 0)
    zero = -99;
  else
    zero = (int) round(min - scalefactor * xmin);
  return (zero);
}

double round(double x)
{
  return(floor((double) x + 0.5));
}

/* plot a vector of doubles */
void plotrvec(FILE *fp, int n, double *y)
{
  int i, xsize;
  double *x=NULL;

  if(x==NULL || xsize < n)
  {
    if(x!=NULL)
      free(x);
    xsize = n;
    x = (double *) mustalloc(xsize * sizeof(double));
    for(i=0; i<xsize; i++)
      x[i] = i;
  }
  charplot(fp, n, x, y);
}

/* plot a vector of integers */
void plotivec(FILE *fp, int n, int *y_in)
{
  int i, vecsize;
  double *x=NULL, *y;

  if(x==NULL || vecsize < n)
  {
    if(x!=NULL)
    {
      free(x);
      free(y);
    }
    vecsize = n;
    x = (double *) mustalloc(vecsize * sizeof(double));
    y = (double *) mustalloc(vecsize * sizeof(double));
    for(i=0; i<vecsize; i++)
      x[i] = i;
  }
  for(i=0; i<vecsize; i++)
    y[i] = y_in[i];
  charplot(fp, n, x, y);
}
