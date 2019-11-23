/*
 * Driver for getf, which lives in identity.c 
 */
#include <stdio.h>
#include "mismatch.h"
#define VECSIZE 15
void main(void)
{
  int  i;
  double f[VECSIZE], *rval;
  double theta[2], tau[1];

  theta[0] = 1.0;
  theta[1] = 10.0;
  tau[0]    = 4.0;

  init_theory(VECSIZE, tau, 1, theta, 2);

  rval = getf(f, VECSIZE, theta[0], theta[1], tau[0]);

  printf("rval=%u f=%u (should be equal)\n", (unsigned) rval,
	 (unsigned) f);
  

  printf("f:");
  for(i=0; i<VECSIZE; i++)
    printf(" %f", f[i]);
  putchar('\n');
  exit(0);
}
