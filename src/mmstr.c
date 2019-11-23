#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "alloc2d.h"
#include "version.h"
#include "readstr.h"
#include "mismatch.h"
#include "bye.h"
#include "header.h"
#include "str_misc.h"
#include "unirand.h"
#define MISSING -9
#define MAXGROUPS 10
#define HEIGHT_RATIO 0.7
#define HSEP 0.2
#define VSEP 0.3
#define FIGCOLS 6
#define ORIGIN 0
#define MEAN 1
/*** prototypes ***/
void usage(void);
int dblcompar(double *x, double *y);
double get_variance(int *f, int max, int about);
void get_moments(double *m2,  double *m4, int *f, int max, int about);
void get_theta0_tau(double *theta0, double *tau, int *f, int max);
double dblquantile(double p, double *vec, int size);

#define NQUANT 7

/*** externals ***/
double qval[NQUANT] = {0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975};
int prfreq = 0;
int verbose = 0;
int which = -1;  /* which group? */
int pictex = 0;
int mismatch = 0;
int doestimate = 1;
int boot = 0;    /* how many bootstrap replicates? */

void usage(void)
{
  fflush(stdout);
  fprintf(stderr,"\nusage: mmstr [options] inputfile");
  fprintf(stderr,"\n  where options may include:");
  fprintf(stderr,"\n  -b<x>  number of bootstrap replicates. Def: %d",
	  boot);
  fprintf(stderr,"\n  -m     print mismatch distributions?");
  fprintf(stderr," Default = %s", YES(mismatch));
  fprintf(stderr,"\n  -p     pictex-format output");
  fprintf(stderr,"\n  -w<g>  which group?");
  fprintf(stderr," Default = %d", which);
  putc('\n', stderr);
  exit(1);
}

/* Compare doubles */
int dblcompar(double *x, double *y)
{

  if(*x < *y)
    return(-1);
  if(*x > *y)
    return(1);
  return(0);
}

/* get a quantile from a sorted array of data */
double dblquantile(double p, double *vec, int size)
{
  int i;

  i = p*size - 1;
  if(i<0)
    i = 0;
  return(vec[i]);
}

/* this algorithm is very slow, very stable */
double get_variance(int *f, int max, int about)
{
  int i;
  double sum = 0.0, m1=0.0, var=0.0, d;

  /* get sum */
  for(i=0; i<=max; i++)
    sum += f[i];

  /* get m1 */
  if(about == ORIGIN)
    m1 = 0.0;
  else
  {
#ifndef NDEBUG
    assert(about == MEAN);
#endif
    for(i=0; i<=max; i++)
      m1 += i * (f[i]/sum);
  }

  /* get variance */
  for(i=0; i<=max; i++)
  {
    d = i - m1;
    var += d*d*(f[i]/(sum-1));
  }
  return(var);
}

/* this algorithm is very slow, very stable */
void get_moments(double *m2,  double *m4, int *f, int max, int about)
{
  int i;
  double sum = 0.0, m1=0.0, d;

  /* get sum */
  for(i=0; i<=max; i++)
    sum += f[i];

  /* get m1 */
  if(about == ORIGIN)
    m1 = 0.0;
  else
  {
#ifndef NDEBUG
    assert(about == MEAN);
#endif
    for(i=0; i<=max; i++)
      m1 += i * (f[i]/sum);
  }

  *m2 = *m4 = 0.0;
  /* get moments 2 and 4 */
  if(sum <= 1)
    return;
  for(i=0; i<=max; i++)
  {
    d = i - m1;
    *m2 += d*d*(f[i]/(sum-1));
    *m4 += d*d*d*d*(f[i]/(sum-1));
  }
}

/* estimate theta0 and tau using 2-parameter method of moments */
void get_theta0_tau(double *theta0, double *tau, int *f, int max)
{
  double m2, m4;

  get_moments(&m2, &m4, f, max, ORIGIN);
  *theta0 = (m4 - m2)/3 - m2*m2 ;
  if(*theta0 <= 0.0)
    *theta0 = 0.0;
  else
    *theta0 = sqrt((double) *theta0);
  *tau = m2 - *theta0;
}

void main(int argc, char **argv)
{
  int gid[MAXGROUPS];  /* gid[i] is identifier of i'th group */
  int smp[MAXGROUPS];  /* smp[i] = number of individuals in group i */
  int ngrps, nloci, maxrepeat;
  int ***frq;          /* frq[i][j][k] = # w/ count k in locus j, group i */
  int **mm=NULL;       /* mismatch distribution */
  int i, j, k;
  double theta0, tau, *btheta0, *btau;
  double vmean, v, theta0mean, taumean;
  FILE *fp = stdin;

  /*  Command line arguments */
  for(i=1; i<argc; ++i)
  {
    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
      case 'b':
	boot = strtol(argv[i]+2, NULL, 10);
	break;
      case 'm':
	mismatch = !mismatch;     /* toggle */
	break;
      case 'p':
	pictex = !pictex;
	break;
      case 'w':
	which = atoi(argv[i]+2);
	break;
      default:
	usage();
      }
    }else
    {
      fp = fopen(argv[i],"r");
      if(fp == NULL)
      {
	fprintf(stderr,"\nCan't open file \"%s\".\n", argv[i]);
	exit(1);
      }
      fprintf(stdout,"#Input file: %s\n", argv[i]);
    }
  }
  if(fp == stdin)
    usage();

  header("mmstr", "(Analysis of STR mismatch distributions)", stdout);
  readstr(fp, gid, smp, MAXGROUPS, &frq, &ngrps, &nloci, &maxrepeat);

  /* allocate arrays for mismatch distributions */
  mm = (int **) alloc2d(nloci+1, maxrepeat+1, sizeof(int));
  if(mm == NULL)
    error("alloc3d couldn't allocate memory");
  for(j=0; j<=nloci; j++)
    for(k=0; k<=maxrepeat; k++)
      mm[j][k] = 0;

  /* which populations are we going to analyze? */
  if(which>=0)
  {
    i = findndx(which, gid, ngrps);
    if(i<0)
      error("group specified with -w was not in data");
    fprintf(stdout,"%c Data for group %d:\n", START_COMMENT, which);
    which = i;
  }else   /* ignore distinctions between groups */
  {
    fprintf(stdout,"%c Data for entire population\n", START_COMMENT);
    /* dump all groups into group 0 */
    for(j=0; j<nloci; j++)
      for(k=0; k<=maxrepeat; k++)
      {
	for(i=1; i<ngrps; i++)
	  frq[0][j][k] += frq[i][j][k];
      }
    which = 0;
  }

  /* fill and print mismatch arrays */
  for(i=0; i<=maxrepeat; i++)
    mm[nloci][i] = 0;
  if(mismatch)
    fprintf(stdout,"%% Mismatch Distributions\n");
  vmean = theta0mean = taumean = 0.0;
  for(k=0; k<nloci; k++)
  {
    get_mismatch(mm[k], frq[which][k], frq[which][k], maxrepeat);
    for(i=0; i<=maxrepeat; i++)
      mm[nloci][i] += mm[k][i];
    vmean += v = get_variance(mm[k], maxrepeat, MEAN);
    get_theta0_tau(&theta0, &tau, mm[k], maxrepeat);
    theta0mean += theta0;
    taumean += tau;
    fprintf(stdout,"\n#Locus %d var=%g theta0=%g tau=%g",
	    k, v, theta0, tau);
    if(mismatch)
    {
      fprintf(stdout,"\nmismatch:");
      for(i=0; i<=maxrepeat; i++)
	fprintf(stdout," %d", mm[k][i]);
    }
  }
  vmean /= nloci;
  theta0mean /= nloci;
  taumean /= nloci;
  if(mismatch)
  {
    fprintf(stdout,"\n%%Aggregate MM dist across loci:\n");
    for(i=0; i<=maxrepeat; i++)
      fprintf(stdout," %d", mm[nloci][i]);
  }
  fprintf(stdout,"\nMeans across %d loci: var=%g theta0=%g tau=%g",
	  nloci, vmean, theta0mean, taumean);

  if(boot)      /* do bootstrap */
  {
    btheta0 = (double *) mustalloc( boot * sizeof(double));
    btau = (double *) mustalloc( boot * sizeof(double));
    for(i=0; i<boot; i++)
    {
      btheta0[i] = btau[i] = 0.0;
      for(j=0; j<nloci; j++)
      {
	k = randint(nloci);  /* random locus */
	get_theta0_tau(&theta0, &tau, mm[k], maxrepeat);
	btheta0[i] += theta0;
	btau[i] += tau;
      }
      btheta0[i] /= nloci;
      btau[i] /= nloci;
    }
    qsort(btheta0, (unsigned) boot, sizeof(double),
	  (int (*)(const void*,const void*)) dblcompar);
    qsort(btau, (unsigned) boot, sizeof(double),
	  (int (*)(const void*,const void*)) dblcompar);
    fprintf(stdout,"\n\nDid %d bootstrap replicates", boot);
    fprintf(stdout,"\n%8s %15s %15s", "Quantile", "theta0", "tau");
    for(i=0; i<NQUANT; i++)
      fprintf(stdout,"\n%8.4f %15.9g %15.9g",
	      qval[i],
	      dblquantile(qval[i], btheta0, boot),
	      dblquantile(qval[i], btau, boot));
  }
  putc('\n', stdout);
  exit(0);
}







