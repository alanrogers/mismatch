#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "eprintf.h"
#include "mytypes.h"
#include "iscoales.h"
#include "mmtheory.h"
/*
 * Exercise functions in mmtheory.  Produces no output if all goes well, 
 * unless the -v (verbose) flag is used.
 */

/* prototypes */
void xgetpoisson(unsigned dim, double tau, unsigned verbose, double *truval);
void xgeteq(unsigned dim, double theta, unsigned verbose);
void xnextf(unsigned dim, double theta, double tau, unsigned verbose);
void xgetf(unsigned dim, double theta0, double theta1, double tau,
	   unsigned verbose);

enum {MAXDIM = 10000};

void main(int argc, char **argv)
{
  int i;
  unsigned verbose = 0;
  double *truval;

  for(i=1; i<argc; i++)
    if(argv[i][0] == '-')
      switch(argv[i][1]) {
      case 'v':
	verbose = !verbose;
	break;
      default:
	fprintf(stderr, "usage: xmmtheory [-v]\n");
	fprintf(stderr, " Options:\n");
	fprintf(stderr, "  -v : verbose output\n");
	exit(2);
      }

  /* allocate truval vector w/ maximal dimension */
  truval = (double *) emalloc(MAXDIM * sizeof(double));
  truval[0] = 1.0;
  for(i=1; i<MAXDIM; i++)
    truval[i] = 0.0;

  /* test getpoisson */
  xgetpoisson(0, 0.0, verbose, truval);
  xgetpoisson(1, 0.0, verbose, truval);
  xgetpoisson(2, 0.0, verbose, truval);
  xgetpoisson(10, 0.0, verbose, truval);
  xgetpoisson(100, 0.0, verbose, truval);

  truval[0] = exp(-2.0);
  for(i=1; i<MAXDIM; i++)
    truval[i] = truval[i-1] * 2.0 / i;
  xgetpoisson(0, 2.0, verbose, truval);
  xgetpoisson(1, 2.0, verbose, truval);
  xgetpoisson(2, 2.0, verbose, truval);
  xgetpoisson(10, 2.0, verbose, truval);
  xgetpoisson(100, 2.0, verbose, truval);

  free(truval);

  /* test geteq */
  xgeteq(0, 0.0, verbose);
  xgeteq(1, 0.0, verbose);
  xgeteq(2, 0.0, verbose);
  xgeteq(10, 0.0, verbose);
  xgeteq(100, 0.0, verbose);

  xgeteq(0, 0.1, verbose);
  xgeteq(1, 0.1, verbose);
  xgeteq(2, 0.1, verbose);
  xgeteq(10, 0.1, verbose);
  xgeteq(100, 0.1, verbose);

  xgeteq(0, 2.0, verbose);
  xgeteq(1, 2.0, verbose);
  xgeteq(2, 2.0, verbose);
  xgeteq(10, 2.0, verbose);
  xgeteq(100, 2.0, verbose);

  /* text nextf */
  xnextf(0, 0.0, 3.0, verbose);
  xnextf(1, 0.0, 3.0, verbose);
  xnextf(2, 0.0, 3.0, verbose);
  xnextf(10, 0.0, 3.0, verbose);
  xnextf(100, 0.0, 3.0, verbose);
  
  xnextf(0, 2.0, 3.0, verbose);
  xnextf(1, 2.0, 3.0, verbose);
  xnextf(2, 2.0, 3.0, verbose);
  xnextf(10, 2.0, 3.0, verbose);
  xnextf(100, 2.0, 3.0, verbose);

  /* text getf */
  xgetf(  0, 1.0, 100.0, 3.0, verbose);
  xgetf(  1, 1.0, 100.0, 3.0, verbose);
  xgetf(  2, 1.0, 100.0, 3.0, verbose);
  xgetf( 10, 1.0, 100.0, 3.0, verbose);
  xgetf(100, 1.0, 100.0, 3.0, verbose);

  xgetf(100,   0.0,   0.0, 3.0, verbose);
  xgetf(100,   0.0, 100.0, 3.0, verbose);
  xgetf(100, 100.0,   0.0, 3.0, verbose);
  xgetf(100, 100.0, 100.0, 3.0, verbose);
  xgetf(100, 0.01,   4.0, 3.0, verbose);
  xgetf(100, 4.0,   0.01, 3.0, verbose);
  xgetf(100, 2.0, 20.0, 0.0, verbose);
  xgetf(100, 2.0, 20.0, 1.0, verbose);
  xgetf(100, 2.0, 20.0, 5.0, verbose);
  xgetf(100, 2.0, 20.0, 10.0, verbose);
  xgetf(100, 2.0, 20.0, 20.0, verbose);
  xgetf(100, 2.0, 20.0, 40.0, verbose);
  xgetf(300, 2.0, 20.0, 80.0, verbose);
  xgetf(800, 2.0, 20.0, 200.0, verbose);
  xgetf(10000, 2.0, 20.0, 1000.0, verbose);
  xgetf(100, 0.000001, 0.000001, 10.0, verbose);
  xgetf(100, 0.000001, 1.0, 10.0, verbose);
  xgetf(100, 1.0, 0.000001, 10.0, verbose);

  exit(0);
}

void xgetpoisson(unsigned dim, double tau, unsigned verbose, double *truval)
{
  double *poi, *poi2, cum, mean, var;
  unsigned i;

  if(dim > MAXDIM)
    eprintf("error in xgetpoisson: dim=%u but MAXDIM=%d", dim, MAXDIM);
  if(verbose)
    printf("xgetpoisson: dim = %u tau = %g\n", dim, tau);
  poi = (double *) emalloc(dim * sizeof(double));
  poi2 = getpoisson(poi, dim, tau);
  if(dim == 0) {
    if(poi2 != NULL)
      eprintf(" getpoisson returned %x.  Should be NULL", (unsigned) poi2);
  }else {
    if(poi2 != poi)
      eprintf(" getpoisson returned %x.  Should be %x",
	      (unsigned) poi2, (unsigned) poi);
  }
  cum = mean = var = 0.0;
  for(i=0; i<dim; i++) {
    assert(poi[i] >= 0.0);
    mean += i * poi[i];
    var += i*i*poi[i];
    cum += poi[i];
    if(verbose && (i <= 10 || i >= dim-11))
      printf(" poi[%3u] = %g  cum = %g\n", i, poi[i], cum);
    if(truval[i] >=0.0 && poi[i] != truval[i])
      printf(
	" discrepancy in xgetpoisson(dim=%u, tau=%g): poi[%u]=%g, not %g\n",
	dim, tau, i, poi[i], truval[i]);
  }
  if(fabs(cum-1.0) < 0.001) {
    var -= mean*mean;
    if(verbose || (fabs(mean-tau) > 0.001 || fabs(var-tau) > 0.01))
      printf(" xgetpoisson(dim=%u,tau=%g): mean = %g var = %g\n",
	     dim, tau, mean, var);
  }
  if(poi != NULL)
    free(poi);
}

void xgeteq(unsigned dim, double theta, unsigned verbose)
{
  double *eq, *eq2, cum, mean, var, evar, truval, tol;
  unsigned i;

  /* theoretical mean = theta. theoretical variance = theta + theta^2 */
  evar = theta + theta*theta;

  if(dim > MAXDIM)
    eprintf("error in xgeteq: dim=%u but MAXDIM=%d", dim, MAXDIM);
  if(verbose)
    printf("xgeteq: dim = %u theta = %g\n", dim, theta);
  eq = (double *) emalloc(dim * sizeof(double));
  eq2 = geteq(eq, dim, theta);
  if(dim == 0) {
    if(eq2 != NULL)
      eprintf(" geteq returned %x.  Should be NULL", (unsigned) eq2);
  }else {
    if(eq2 != eq)
      eprintf(" geteq returned %x.  Should be %x",
	      (unsigned) eq2, (unsigned) eq);
  }
  cum = mean = var = 0.0;
  for(i=0; i<dim; i++) {
    assert(eq[i] >= 0.0);
    mean += i * eq[i];
    var += i*i*eq[i];
    cum += eq[i];
    if(verbose && (i <= 10 || i >= dim-11))
      printf(" eq[%3u] = %g  cum = %g\n", i, eq[i], cum);
    truval = pow(theta, (double) i) / pow(1.0+theta, i+1.0);
    if(truval > 1.0)
      tol = truval * DBL_EPSILON;
    else
      tol = DBL_EPSILON;
    if(fabs(eq[i]-truval) > tol) {
      printf(" discrepancy in xgeteq(dim=%u, theta=%g): eq[%u]=%g, not %g.",
	dim, theta, i, eq[i], truval);
      printf(" diff = %g\n", eq[i] - truval);
    }

  }
  if(fabs(cum-1.0) < 0.001) {
    var -= mean*mean;
    if(verbose || (fabs(mean-theta) > 0.001))
      printf(" xgeteq(dim=%u,theta=%g): mean = %g (should be %g)\n",
	     dim, theta, mean, theta);
    if(verbose || (fabs(var-evar) > 0.01)) 
      printf(" xgeteq(dim=%u,theta=%g): var = %g (should be %g)\n",
	     dim, theta, var, evar);
  }
  if(eq != NULL)
    free(eq);
}

void xnextf(unsigned dim, double theta, double tau, unsigned verbose)
{
  double *f0, *f, *g, relerr;
  unsigned i;

#if 0
  /* theoretical moments about zero */
  A = (theta0 - theta1)*exp(-tau/theta1);
  mean = A + theta1;
  u2 = theta1 + 2.0*theta1*theta1 
    + A*(1.0 + 2.0*tau + 2.0*(theta0 + theta1));
  u3 = theta1 + 6.0*theta1*theta1*(1.0 + theta1)
    + A*(1.0 + 3.0*tau*tau + 6.0*((tau + theta0)*(1.0 + theta0 + theta1)
			    +theta1*(1.0 + theta1)));
#endif

  if(dim > MAXDIM)
    eprintf("error in xnextf: dim=%u but MAXDIM=%d", dim, MAXDIM);
  if(verbose)
    printf("xnextf: dim = %u theta = %g tau=%g\n", dim, theta, tau);
  f0 = (double *) emalloc(dim * sizeof(double));
  f  = (double *) emalloc(dim * sizeof(double));

  /* f0 is equilibrium (theta) */
  g = geteq(f0, dim, theta);
  if(dim == 0) {
    if(g != NULL)
      eprintf(" geteq returned %x.  Should be NULL", (unsigned) g);
  }else {
    if(g != f0)
      eprintf(" geteq returned %x.  Should be %x",
	      (unsigned) g, (unsigned) f0);
  }

  /* new f has same theta, so distribution should be unchanged */
  g = nextf(f, f0, dim, theta, tau);
  if(dim == 0) {
    if(g != NULL)
      eprintf(" nextf returned %x.  Should be NULL", (unsigned) g);
  }else {
    if(g != f)
      eprintf(" geteq returned %x.  Should be %x",
	      (unsigned) g, (unsigned) f);
  }

  for(i=0; i<dim; i++) {
    assert(f[i] >= 0.0);
    if(verbose && (i <= 10 || i >= dim-11))
      printf(" f0[%3u] = %g  f[%3u] = %g\n", i, f0[i], i, f[i]);
    relerr = fabs(f[i] - f0[i])/f0[i];
    if(relerr > DBL_EPSILON) {
      printf(" discrepancy in xnextf(dim=%u, theta=%g): f[%u]=%g, not %g.",
	dim, theta, i, f[i], f0[i]);
      printf(" diff = %g relerr = %g\n", f[i] - f0[i], relerr);
    }
  }
  if(f0 != NULL)
    free(f0);
  if(f != NULL)
    free(f);
}

void xgetf(unsigned dim, double theta0, double theta1, double tau,
	   unsigned verbose)
{
  double *f, *g, relerr, A, u1, u2, u3, Eu1, Eu2, Eu3, cum=0.0, tol;
  unsigned i, identified=0;

  if(verbose) {
    printf("xgetf: (dim,theta0,theta1,tau)=(%u,%g,%g,%g)\n",
	   dim,theta0,theta1,tau);
    identified=1;
  }

  if(dim > MAXDIM)
    eprintf("error in xgetf: dim=%u but MAXDIM=%d", dim, MAXDIM);
  f  = (double *) emalloc(dim * sizeof(double));

  /* mismatch distribution */
  g = getf(f, dim, theta0, theta1, tau);
  if(dim == 0) {
    if(g != NULL)
      eprintf(" getf returned %x.  Should be NULL", (unsigned) g);
  }else {
    if(g != f)
      eprintf(" getf returned %x.  Should be %x",
	      (unsigned) g, (unsigned) f);
  }

  /* calculate moments about zero */
  u1 = u2 = u3 = 0.0;
  for(i=0; i<dim; i++) {
    cum += f[i];
    u1 += i * f[i];
    u2 += i*i * f[i];
    u3 += i*i*i * f[i];
  }

  if(cum < 0.99) {  /* early return */
    if(verbose)
      printf(" early return (cum=%g)\n", cum);
    return;
  }

  /* theoretical moments about zero */
  A = (theta0 - theta1)*exp(-tau/theta1);
  Eu1 = A + theta1;
  Eu2 = theta1 + 2.0*theta1*theta1 
    + A*(1.0 + 2.0*tau + 2.0*(theta0 + theta1));
  Eu3 = theta1 + 6.0*theta1*theta1*(1.0 + theta1)
    + A*(1.0 + 3.0*tau*tau + 6.0*((tau + theta0)*(1.0 + theta0 + theta1)
			    +theta1*(1.0 + theta1)));
  tol = FLT_EPSILON;

  if(verbose) {
    printf(" u1=%g Eu1=%g\n", u1, Eu1);
    printf(" u2=%g Eu2=%g\n", u2, Eu2);
    printf(" u3=%g Eu3=%g\n", u3, Eu3);
    printf(" tol = %g\n", tol);
  }

  /* compare with theoretical values */
  relerr = fabs(u1 - Eu1)/Eu1;
  if(relerr > tol) {
    if(!identified) {
      printf("xgetf: (dim,theta0,theta1,tau)=(%u,%g,%g,%g)\n",
	     dim,theta0,theta1,tau);
      identified = 1;
    }
    printf(" u1=%g, not %g.", u1, Eu1);
    printf(" diff = %g relerr = %g\n", u1 - Eu1, relerr);
  }
  relerr = fabs(u2 - Eu2)/Eu2;
  if(relerr > tol) {
    if(!identified) {
      printf("xgetf: (dim,theta0,theta1,tau)=(%u,%g,%g,%g)\n",
	     dim,theta0,theta1,tau);
      identified = 1;
    }
    printf(" u2=%g, not %g.", u2, Eu2);
    printf(" diff = %g relerr = %g\n", u2 - Eu2, relerr);
  }
  relerr = fabs(u3 - Eu3)/Eu3;
  if(relerr > tol) {
    if(!identified) {
      printf("xgetf: (dim,theta0,theta1,tau)=(%u,%g,%g,%g)\n",
	     dim,theta0,theta1,tau);
      identified = 1;
    }
    printf(" u3=%g, not %g.", u3, Eu3);
    printf(" diff = %g relerr = %g\n", u3 - Eu3, relerr);
  }
  if(f != NULL)
    free(f);
}

