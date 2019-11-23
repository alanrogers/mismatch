#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "bye.h"
#define TOGGLE(i) ((i) == 0)
#define YES(i) ((i) ? "YES" : "NO")
/****************************************************************
Normalize a distribution so that it sums to 1, or invert this
process to obtain a distribution of counts.
****************************************************************/
/** prototypes **/
void usage(void);
char *mustalloc(unsigned bytes);
void option(char *flag, char *msg, char *def);

/*** external variables **/
int sampsize = 0;
int twocol = 0;

/** describe an option w/i usage() function **/
void option(char *flag, char *msg, char *def)
{
  fprintf(stderr,"\n %-5s %s Default=%s", flag, msg, def);
}
/** allocate memory using malloc; abort on failure **/
char *mustalloc(unsigned bytes)
{
  char *p;

  p = malloc(bytes);
  if(p==NULL)
  {
    fflush(stdout);
    fprintf(stderr,"\nmustalloc: Can't allocate %d bytes of memory.\n",
	    bytes);
    exit(1);
  }
  return(p);
}
int main(int argc, char **argv)
{
  double sum, *f, x, y;
  int i, dim;
  FILE *fp = NULL;

  /* process command line arguments */
  for(i=1; i<argc; i++)
  {
    if(argv[i][0] == '-')
      switch(argv[i][1])
      {
      case 'r':
	sampsize = atoi(argv[i]+2);
	break;
      case '2':
	twocol = TOGGLE(twocol);
	break;
      default:
	usage();
      }
    else
    {
      fp = fopen(argv[i], "r");
      if(fp == NULL)
      {
	fprintf(stderr,"\nCan't read file %s\n", argv[i]);
	exit(1);
      }
    }
  }
  if(fp == NULL)
    usage();

  /* count observations */
  for(dim=0;
      fscanf(fp, "%lf", &sum)==1;
      dim++)
    ;
  if(dim < 2)
  {
    fprintf(stderr,"\nError: only %d item(s) of data\n", dim);
    exit(1);
  }
  /* allocate array */
  f = (double *) mustalloc(dim * sizeof(double));

  /* read data, accumulate sum */
  rewind(fp);
  sum = 0.0;
  for(i=0; i<dim; i++)
  {
    if(fscanf(fp,"%lf", f+i) != 1)
    {
      fprintf(stderr,"\nBad input\n");
      exit(1);
    }
    if(!sampsize)
      sum += f[i];
  }
  if(sampsize)
    sum = sampsize*(sampsize-1)/2;
  else  /*  report n */
  {
    x = 0.5*(1.0 + sqrt(1.0+8.0*sum));
    y = floor(x + 0.5);  /* round to nearest int */
    if(fabs(x-y) < 0.001)
      printf("\nN = %.0f", y);
    else
      printf("\nN = %f...[Warning: this should be an integer]",
	     x);
  }
  for(i=0; i<dim; i++)
  {
    if(sampsize)
    {
      y = f[i]*sum;
      x = floor(y+0.5);
      if(fabs(x-y) > 0.5)
      {
	fprintf(stderr,"\nError: sampsize can't be %d\n", sampsize);
	exit(1);
      }
      if(twocol)
	printf("\n%d %.0f", i, x);
      else
	printf("\n%.0f", x);
    }else
    {
      x = f[i]/sum;
      if(twocol)
	printf("\n%d %f", i, x);
      else
	printf("\n%f", x);
    }
  }
  putc('\n', stdout);
  return 0;
}
void usage(void)
{
  fprintf(stderr,"\nusage: normdist [options] inputfile");
  fprintf(stderr,"\n  where options may include:");
  option("-r<sampsize>", "Invert the usual operation?", YES(sampsize));
  option("-2", "Two-column output?", YES(twocol));
  putc('\n', stderr);
  exit(1);
}
  
