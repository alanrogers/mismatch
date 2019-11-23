#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "memstack.h"
#include "mismatch.h"
#include "getcic.h"
#include "bye.h"
#define BUFFSIZE 200
#define VECSIZE 200
#define WORDSIZE 100
/** Definitions relating to keywords appear only in this file **/
enum keyword getkeyword(FILE *fp);


/****************************************************************
Input routine.  Returns allocated length of histogram.
****************************************************************/
int    input(SIMULATION *s, double *from, double *to, double *by, int msize,
	     unsigned *eflags, FILE *fp, FILE *ofp)
{
  int  i, ival;
  ASSIGNMENT *a;


  *eflags = 0;
  while((a=getassignment(fp)) != NULL)
  {
    switch(a->lhs)
    {
    case InputFile:
      s->fname = dupstring(a->val[0]);
      break;
    case Sampsize:
      s->smpsiz = atoi(a->val[0]);
      break;
    case NSites:
      s->nsites = atoi(a->val[0]);
      break;
    case Histogram:
      if(msize <= 0)   /* histogram size set by input data */
	msize = a->length;
      /* otherwise, histogram size is set by parameter msize */
      s->h = (int *) mustalloc(msize*sizeof(int));
      s->h[msize-1] = 0;
      fprintf(ofp,
      "\n%c Using histogram of length %d. Input histogram had length %d.",
	   START_COMMENT, msize, a->length);
      for(i=0; i<a->length; i++)
      {
	ival = atoi(a->val[i]);  /* round to nearest int */
	/* ensure that input value was an int */
	if(fabs(((double) ival) - atof(a->val[i])) > 0.000001)
	{
	  fflush(stdout);
	  fprintf(stderr,"\nError in input: Non-integer value in histogram");
	  fprintf(stderr,"\n  Illegal value: h[%d] = %s\n",
		  i, a->val[i]);
	  exit(1);
	}
	if(i >= msize)           /* last entry of h is sum of large values */
	  s->h[msize-1] += ival;
	else
	  s->h[i] = ival;
      }
      for(i=a->length; i<msize; i++)  /* zero unassigned portion of h */
	s->h[i] = 0;
      break;
    case Cumulants:
      assert(a->length == 3);
      for(i=0; i<3; i++)
	s->e.cumulant[i] = atof(a->val[i]);
      break;
    case Labels:
      for(i=0; i < a->length; i++)
      {
	if(!strcmp(lowercase(a->val[i]), "theta0"))
	  *eflags |= THETA0;
	if(!strcmp(lowercase(a->val[i]), "theta1"))
	  *eflags |= THETA1;
	if(!strcmp(lowercase(a->val[i]), "tau"))
	  *eflags |= TAU;
	if(!strcmp(lowercase(a->val[i]), "mse"))
	  *eflags |= MSE;
	if(!strcmp(lowercase(a->val[i]), "mae"))
	  *eflags |= MAE;
	if(!strcmp(lowercase(a->val[i]), "roughness"))
	  *eflags |= ROUGHNESS;
	if(!strcmp(lowercase(a->val[i]), "seg"))
	  *eflags |= SEG;
      }
      break;
    case Estimates:
      assert(a->length == nestimates);
      i=0;
      if(*eflags & THETA0)
	s->e.theta0 = atof(a->val[i++]);
      if(*eflags & THETA1)
	s->e.theta1 = atof(a->val[i++]);
      if(*eflags & TAU)
	s->e.tau = atof(a->val[i++]);
      if(*eflags & MSE)
	s->e.mse = atof(a->val[i++]);
      if(*eflags & MAE)
	s->e.mae = atof(a->val[i++]);
      if(*eflags & ROUGHNESS)
	s->e.roughness = atof(a->val[i++]);
      if(*eflags & SEG)
	s->e.seg = atof(a->val[i++]);
      assert(i == nestimates);
      break;
    case Seg:
      s->e.seg = atof(a->val[0]);
      break;
    case RangeLog10Theta0:
      if(from==NULL || to==NULL || by==NULL)
	break;
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[0] = atof(a->val[0]);
      to[0]   = atof(a->val[1]);
      by[0]   = atof(a->val[2]);
      break;
    case RangeGrowth:
      if(from==NULL || to==NULL || by==NULL)
	break;
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[1] = atof(a->val[0]);
      to[1]   = atof(a->val[1]);
      by[1]   = atof(a->val[2]);
      break;
    case RangeTau:
      if(from==NULL || to==NULL || by==NULL)
	break;
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[2] = atof(a->val[0]);
      to[2]   = atof(a->val[1]);
      by[2]   = atof(a->val[2]);
      break;
    case Eof:
      break;
    case Badval:
    default:
      error("input(): Unrecognized keyword in input");
    }
    freeassignment(a);
  }
  return (msize);
}
/****************************************************************
Write a floating point variable's name and value to output
****************************************************************/
void fputval(char *name, int len, double *vec, FILE *fp, int comment)
{
  putc('\n', fp);
  if(comment)
    putc(START_COMMENT, fp);
  fprintf(fp,"%s =", name);
  while(len-- > 0)
    fprintf(fp, " %9g", *vec++);
  putc(';', fp);
}
/****************************************************************
Write an int variable's name and value to output
****************************************************************/
void iputval(char *name, int len, int *vec, FILE *fp, int comment)
{
  putc('\n', fp);
  if(comment)
    putc(START_COMMENT, fp);
  fprintf(fp,"%s =", name);
  while(len-- > 0)
    fprintf(fp, " %d", *vec++);
  putc(';', fp);
}
/****************************************************************
Write a (string) variable's name and value to output
****************************************************************/
void sputval(char *name, int len, char **vec, FILE *fp, int comment)
{
  putc('\n', fp);
  if(comment)
    putc(START_COMMENT, fp);
  fprintf(fp,"%s =", name);
  while(len-- > 0)
    fprintf(fp, " %s", *vec++);
  putc(';', fp);
}
/****************************************************************
Read an assignment  from input file line of form:
  name = val1 val2 val3 ... valK ;
Return a structure containing left hand side of assignment and a vector
containing the right hand side.  Return NULL on end of file.  
****************************************************************/
ASSIGNMENT *getassignment(FILE *fp)
{
  ASSIGNMENT *a;
  char buff[BUFFSIZE];
  char *vec[VECSIZE];
  int vecsize, i;
  enum keyword word;

  /** Allocate and initialize assignment **/
  a = (ASSIGNMENT *) mustalloc(sizeof(ASSIGNMENT));
  a->lhs = Badval;
  a->val = NULL;
  a->length = 0;

  a->lhs = getkeyword(fp);  /* read lhs of equation */
  if(a->lhs == Eof)  /* end of file */
  {                  /* This is the only place where EOF is legal */
    free(a);
    return(NULL);
  }

  word = getkeyword(fp);      /* read '=' */
  if(word != Equals)       /* make sure an "=" was read */
    error("Missing equal sign in input data");

  /* Right hand side is a vector of unknown length, terminated by ';' */
  for(vecsize=0; ; vecsize++)
  {
    if(getstring(buff,BUFFSIZE,fp) == NULL)  /* end of file */
      error("Unexpected EOF on input.  Did you forget a semicolon?");
    if(strcmp(buff,";")==0)                  /* semicolon */
      break;
    if(vecsize >= VECSIZE)        /* check for oversized input vector */
    {
      fprintf(stderr,"\nERROR: Input vector is too long.  Max=%d\n",
	      VECSIZE);
      exit(1);
    }
    vec[vecsize] = dupstring(buff);
  }
  if(vecsize == 0)
    error("RHS of assignment is empty in input data.");
  a->length = vecsize;
  a->val = (char **) mustalloc(vecsize*sizeof(char *));
  for(i=0; i<vecsize; i++)
    a->val[i] = vec[i];
  return(a);
}

/****************************************************************
Read a keyword, return numeric code.  This code is very slow, but that 
isn't worth fixing because it is used only briefly at startup.
****************************************************************/
enum keyword getkeyword(FILE *fp)
{
  char buff[BUFFSIZE];
  int i;

  if(getstring(buff, BUFFSIZE, fp) == NULL)
    return(Eof);

  for(i=0; buff[i] != '\0'; i++)
    buff[i] = tolower(buff[i]);  /* ignore case */

  if(strcmp(buff, ";") == 0)
    return(Semicolon);
  if(strcmp(buff, "=") == 0)
    return(Equals);
  if(strcmp(buff, "test") == 0)
    return(Test);
  if(strcmp(buff, "inputfile") == 0)
    return(InputFile);
  if(strcmp(buff, "nsequences") == 0)
    return(Sampsize);
  if(strcmp(buff, "nsites") == 0)
    return(NSites);
  if(strcmp(buff, "mismatch") == 0)
    return(Histogram);
  if(strcmp(buff, "cumulants") == 0)
    return(Cumulants);
  if(strcmp(buff, "labels") == 0)
    return(Labels);
  if(strcmp(buff, "estimates") == 0)
    return(Estimates);
  if(strcmp(buff, "segregating_sites") == 0)
     return(Seg);
  if(strcmp(buff, "rangelog10theta0") == 0)
    return(RangeLog10Theta0);
  if(strcmp(buff, "rangegrowth") == 0)
    return(RangeGrowth);
  if(strcmp(buff, "rangetau") == 0)
    return(RangeTau);
  if(strcmp(buff, "growth") == 0)
    return(Growth);

  /* If no matching keyword is found, return Badval */
  fprintf(stderr,"\nWarning: \"%s\" is not a keyword.",
	  buff);
  return(Badval);
}

/****************************************************************
Read a character string from a file.  Strings are delimited by white space. 
In addition, "=" and ";" are put into separate strings.
****************************************************************/
char *getstring(char *str, int maxlen, FILE *fp)
{
  int c, i;
  
  do{
    c = getcic(fp);
  }while(isspace(c) && c!=EOF);  /* skip initial whitespace */
  if(c==EOF)
    return(NULL);

  if(c==';' || c=='=')        /* ";" and "=" form 1-char strings */
  {
    str[0] = c;
    str[1] = '\0';
    return(str);
  }
  i=0;
  do{                        /* otherwise, accumulate chars into str */
    str[i++] = c;
    c = getcic(fp);
  }while(i<maxlen && instring(c));
  if(i == maxlen && instring(c))
    exit(1);
  ungetc(c, fp);
  str[i] = '\0';
  return(str);
}

/* Read a float from a file, ignoring comments */
int getdouble(double *x, FILE *fp)
{
  char buff[100], *b;

  b = getstring(buff, 100, fp);
  if(b==NULL)
    return(EOF);
  *x = atof(buff);
  return(EOF+1);
}

/****************************************************************
instring(c)==1 if c is a "normal" character, 0 if it is a
space, a '=', a ';', or EOF.
****************************************************************/
int instring(int c)
{
  if(isspace(c) || c=='=' || c==';' || c==EOF)
    return(0);
  return(1);
}
    
void freeassignment(ASSIGNMENT *a)
{
  int i;

  if(a==NULL)
    return;
  for(i=0; i<a->length; i++)
    free(a->val[i]);
  if(a->length > 0)
    free(a->val);
  free(a);
}
  
