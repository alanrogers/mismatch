#include "mytypes.h"
#include "bye.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/****************************************************************
This file contains short functions of that print error messages
and abort.  They are in a separate file so that they can be easily
shared among different directories.
****************************************************************/

/** open a file, abort if unsuccessful **/
FILE *mustopen(char *name, char *mode)
{
  FILE *fp;

  fp = fopen(name, mode);
  if(fp==NULL)
  {
    fprintf(stderr,"\nCan't open file \"%s\" with mode \"%s\"\n",
	    name, mode);
    exit(1);
  }
  return(fp);
}

/** print an error message and quit **/
void error(char *s)
{
  fflush(stdout);
  fprintf(stderr,"\nERROR: %s\n", s);
  exit(1);
}

/** allocate memory using malloc; abort on failure **/
char *mustalloc(unsigned bytes)
{
  char *p;

  p = (char *) malloc(bytes);
  if(p==NULL)
  {
    fflush(stdout);
    fprintf(stderr,"\nmustalloc: Can't allocate %d bytes of memory.\n",
	    bytes);
    exit(1);
  }
  return(p);
}

/****************************************************************
Allocate a vector of floats whose index runs from low to high.  On
successful completion, return pointer to new structure of type
FloatVec.  On return, abort.
****************************************************************/ 
FloatVec *newFloatVec(int lo, int hi)
{
  FloatVec *v;
  int len;

  v = (FloatVec *) malloc( sizeof(FloatVec) );
  if(v == NULL)
    error("newFloatVec: out of memory");
  len = hi - lo + 1;
  v->b = v->f = (float *) malloc(len * sizeof(float));
  if(v->f == NULL)
    error("newFloatVec: out of memory");
  v->f -= lo;  /* now f[lo] is 1st entry in vector */
  v->chopLo = v->lo = lo;
  v->chopHi = v->hi = hi;
  return(v);
}

/****************************************************************
Allocate a vector of doubles's whose index runs from low to high.  On
successful completion, return pointer to new structure of type
DblVec.  On return, abort.
****************************************************************/ 
DblVec *newDblVec(int lo, int hi)
{
  DblVec *v;
  int len;

  v = (DblVec *) malloc( sizeof(DblVec) );
  if(v == NULL)
    error("newDblVec: out of memory");
  len = hi - lo + 1;
  v->b = v->f = (double *) malloc(len * sizeof(double));
  if(v->f == NULL)
    error("newDblVec: out of memory");
  v->f -= lo;  /* now f[lo] is 1st entry in vector */
  v->chopLo = v->lo = lo;
  v->chopHi = v->hi = hi;
  return(v);
}

/****************************************************************
Allocate a vector of ints whose index runs from low to high.  On
successful completion, return pointer to new structure of type
Intvec.  On return, abort.
****************************************************************/ 
IntVec *newIntVec(int lo, int hi)
{
  IntVec *v;
  int len;

  v = (IntVec *) malloc( sizeof(IntVec) );
  if(v == NULL)
    error("newIntVec: out of memory");
  len = hi - lo + 1;
  v->b = v->f = (int *) malloc(len * sizeof(int));
  if(v->f == NULL)
    error("newIntVec: out of memory");
  v->f -= lo;  /* now f[lo] is 1st entry in vector */
  v->chopLo = v->lo = lo;
  v->chopHi = v->hi = hi;
  return(v);
}

/* initialize FloatVec to zero */
void initFloatVec(FloatVec *v)
{
    v->chopLo = v->lo;
    v->chopHi = v->hi;
    memset(v->b, 0, sizeof(float) * (v->hi - v->lo + 1));
}

/* initialize DboVec to zero */
void initDblVec(DblVec *v)
{
    v->chopLo = v->lo;
    v->chopHi = v->hi;
    memset(v->b, 0, sizeof(double) * (v->hi - v->lo + 1));
}

/* initialize IntVec to zero */
void initIntVec(IntVec *v)
{
    v->chopLo = v->lo;
    v->chopHi = v->hi;
    memset(v->b, 0, sizeof(int) * (v->hi - v->lo + 1));
}


