#include <stdio.h>
#include <math.h>
#include "mismatch.h"
#define VECLEN 30
#define MEAN 10
/*
 * Test function get_mse.
 */
void main(void)
{
  double f[VECLEN], sum=0.0, mse;
  int h[VECLEN];
  int i, d=1;

  for(i=0; i<VECLEN; i++)
    sum += f[i] = exp(-(i-MEAN)*(i-MEAN)/(2.0*5.0));

  for(i=0; i<VECLEN; i++)
    f[i] /= sum;

  for(i=0; i<VECLEN; i++)
  {
    h[i] = 100*f[i] + d;
    d = -d;
  }

  mse = get_mse(f, h, VECLEN);
  printf("mse=%g log10(mse)=%g\n", mse, log10(mse));
  exit(0);
}
