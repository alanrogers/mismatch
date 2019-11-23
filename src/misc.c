#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "memstack.h"
#include "mismatch.h"
#include "bye.h"

/** echo the command line **/
void  echo_cmdline(int argc, char **argv)
{
  int i;
  
  fputs("\n% Cmd line:", stdout);
  for(i=0; i<argc; i++)
    printf(" %s", argv[i]);
}

/** duplicate a string **/
char *dupstring(char *s1)
{
  char *s2;
  int i, len;

  len = strlen(s1);
  if(len==0)
    return(NULL);

  s2 = (char *) mustalloc((len+1)*sizeof(char));
  for(i=0; i<len; i++)
    s2[i] = s1[i];
  s2[len] = '\0';
  return(s2);
}
  
/* get a quantile from a sorted array of data */
double quantile(double p, double *vec, int size)
{
  int i;

  i = p*size - 1;
  if(i<0)
    i = 0;
  return(vec[i]);
}

/* Compare doubles */
int compar(double *x, double *y)
{

  if(*x < *y)
    return(-1);
  if(*x > *y)
    return(1);
  return(0);
}

void prmat(double **m, int rows, int cols, FILE *fp)
{
  int i, j;

  for(i=0; i<rows; i++)
  {
    fprintf(fp,"\n%c ", START_COMMENT);
    for(j=0; j<cols; j++)
      fprintf(fp," %11.8f", m[i][j]);
  }
}
/** describe an option w/i usage() function **/
void option(char *flag, char *msg, char *def)
{
  if(def != NULL)
    fprintf(stderr,"\n %-5s %s Default=%s", flag, msg, def);
  else
    fprintf(stderr,"\n %-5s %s", flag, msg);
}

void complain(char *s, char *e)
{
  fflush(stdout);
  fprintf(stderr,"Inconsistent \"define\"s in mismatch.h.  ");
  fprintf(stderr,"Must set %s=1 when %s=1.\n",
	  e, s);
  exit(1);
}

/****************************************************************  
Get 1st 3 cumulants  from histogram.  On return,
m[0] = mean (or 1st cumulant),
m[1] = average of (i - m[0])^2 (the variance or 2nd cumulant)
m[2] = average of (i - m[0])^3 (the 3rd cumulant)
****************************************************************/
void getcumulants(int *h, int msize, double *m)
{
  int i;
  double x, y, s=0.0;

  m[0] = m[1] = m[2] = 0.0;
  for(i=0; i<msize; i++)
  {
    m[0] += i*h[i];
    s += h[i];
  }
  s = 1.0/s;
  m[0] *= s;  /* divide by sum of h */
  for(i=0; i<msize; i++)
  {
    x = i-m[0];  /* diff btw i and mean */
    y = x*x;     /* squared diff */
    m[1] += y*h[i];
    y *= x;      /* cubed diff */
    m[2] += y*h[i];
  }
  m[1] *= s;  /* divide by sum of h */
  m[2] *= s;  /* divide by sum of h */
}

/* find length of vector excluding terminal zeroes */
int veclength(int *h, int dim)
{
  while(h[dim-1]==0)
    --dim;
  return(dim);
}

/* print a vector of doubles in pictex format, preceded by a label */
void pictex_rvec(FILE *fp, int len, double *v, char *lbl)
{
  int i;
  
  fprintf(fp,"\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%s", lbl);
  fprintf(fp,"\n\\plot");
  for(i=0; i<len; i++)
  {
    fprintf(fp," %d %f", i, v[i]);
    if((i+1)%5 == 0 && i < len-1)
      fprintf(fp,"\n     ");
  }
  fprintf(fp,"\n/");
}

/* print a message surrounded by comment characters */
void bold_comment(char *msg, FILE *fp, int comment_char)
{
  int i, linewidth=15;
  putc('\n', fp);
  for(i=0; i<linewidth;i++)
    putc(comment_char, fp);
  fputs(msg,fp);
  for(i=0; i<linewidth; i++)
    putc(comment_char, fp);
}

/* print estimates from a simulation structure */
#if 0
void print_estimates(SIMULATION *obs, FILE *fp, unsigned eflags, int docomment)
{
  int c;

  /* decide whether to prepend a comment symbol */
  if(docomment==PREPEND_COMMENT)
    c = START_COMMENT;
  else
    c = '\0';

  if(eflags & THETA0)
    fprintf(fp,"\n%c%20s = %g ;", c, "theta0_estimate", obs->e.theta0);

  if(eflags & THETA1)
    fprintf(fp,"\n%c%20s = %g ;", c, "theta1_estimate", obs->e.theta1);

  if(eflags & TAU)
    fprintf(fp,"\n%c%20s = %g ;", c, "tau_estimate", obs->e.tau);

  if(eflags & MSE)
    fprintf(fp,"\n%c%20s = %g ;", c, "MSE", obs->e.mse);

  if(eflags & MAE)
    fprintf(fp,"\n%c%20s = %g ;", c, "MAE", obs->e.mae);

  if(eflags & ROUGHNESS)
    fprintf(fp,"\n%c%20s = %g ;", c, "roughness_estimate", obs->e.roughness);

  if(eflags & SEG)
    fprintf(fp,"\n%c%20s = %d ;", c, "Segregating_sites", obs->e.seg);
}
#else
void print_estimates(SIMULATION *obs, FILE *fp, unsigned eflags, int docomment)
{
  if(docomment==PREPEND_COMMENT)
    fprintf(fp,"\n%c%9s = ", START_COMMENT, "Estimates");
  else
    fprintf(fp,"\n%9s = ", "Estimates");

  if(eflags & THETA0)
    fprintf(fp," %10g", obs->e.theta0);

  if(eflags & THETA1)
    fprintf(fp," %10g", obs->e.theta1);

  if(eflags & TAU)
    fprintf(fp," %10g", obs->e.tau);

  if(eflags & MSE)
    fprintf(fp," %10g", obs->e.mse);

  if(eflags & MAE)
    fprintf(fp," %10g", obs->e.mae);

  if(eflags & ROUGHNESS)
    fprintf(fp," %10g", obs->e.roughness);

  if(eflags & SEG)
    fprintf(fp," %5d", obs->e.seg);
  fprintf(fp, " ;");
}
#endif
void print_labels(FILE *fp, unsigned eflags, int docomment)
{
  if(docomment==PREPEND_COMMENT)
    fprintf(fp,"\n%c%9s = ", START_COMMENT, "Labels");
  else
    fprintf(fp,"\n%9s = ", "Labels");

  if(eflags & THETA0)
    fprintf(fp," %10s", "theta0");

  if(eflags & THETA1)
    fprintf(fp," %10s", "theta1");

  if(eflags & TAU)
    fprintf(fp," %10s", "tau");

  if(eflags & MSE)
    fprintf(fp," %10s", "MSE");

  if(eflags & MAE)
    fprintf(fp," %10s", "MAE");

  if(eflags & ROUGHNESS)
    fprintf(fp," %10s", "roughness");

  if(eflags & SEG)
    fprintf(fp," %5s", "Seg");
  fprintf(fp," ;");
}

/* count the number of bits turned on in an unsigned int */
unsigned countbits(unsigned u)
{
  int nbits = 0;
  
  while(u > 0)
  {
    if(u & 1)
      nbits++;
    u >>= 1;
  }
  return(nbits);
}

/* convert string to lower case */
char *lowercase(char *s)
{
  char *s2;

  for(s2 = s; *s2 != '\0'; s2++)
    *s2 = tolower(*s2);
  return(s);
}
