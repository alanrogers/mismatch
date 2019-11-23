#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "alloc2d.h"
#include "readstr.h"
#define MISSING -9
#define MAXGROUPS 10

/* find index of x in vector v. return -1 on failure */
int findndx(int x, int *v, int len)
{
  int i;

  for(i=0; i<len; i++)
    if(x == v[i])
      return(i);
  return(-1);
}


void readstr(
  FILE *fp,       /* pointer to input stream */
  int *gid,       /* gid[i] is identifier of i'th group */
  int *smp,       /* smp[i] = number of individuals in group i */
  int maxgroups,  /* dimension of arrays smp[] and gid[] */
  int ****frq,     /* (*frq)[i][j][k] = # w/ count k in locus j, group i */
  int *ngrps,      /* *ngrps is the number of groups */
  int *nloci,     /* *nloci is the number of loci */
  int *maxrepeat  /* *maxrepeat is the max repeat score */
 )
{
  int i, j, k, count, grp;

  for(i=0; i<MAXGROUPS; i++)   /* initialize group sample sizes */
    smp[i] = 0;

  /* 1st pass through data: count groups, find maxrepeats */
  fscanf(fp, "%d", nloci);
  if( *nloci <= 0 )
  {
    fflush(stdout);
    fprintf(stderr,"\nBad input: nloci=%d\n", *nloci);
    exit(1);
  }
  *ngrps = *maxrepeat = 0;
  while( fscanf(fp, "%d", &grp) != EOF)
  {
    if( (k = findndx(grp, gid, *ngrps)) >= 0)
      smp[k] += 1;
    else
    {
      gid[*ngrps] = grp;
      smp[*ngrps] = 1;
      *ngrps +=1;
    }
    for(i=0; i < *nloci; i++)  
    {
      if(fscanf(fp, "%d", &j) == EOF)
      {
	fprintf(stderr,"\nUnexpected EOF.\n");
	exit(1);
      }
      if(j > *maxrepeat)
	*maxrepeat = j;
    }
  }

  /* allocate frequency array */
  *frq = (int ***) alloc3d(*ngrps, *nloci, *maxrepeat+1, sizeof(int));
  if(*frq == NULL)
  {
    fprintf(stderr,"\nERROR: out of memory in readstr.\n");
    exit(1);
  }
  for(i=0; i < *ngrps; i++)
    for(j=0; j < *nloci; j++)
      for(k=0; k <= *maxrepeat; k++)
	(*frq)[i][j][k] = 0;

  /* 2nd pass through data: fill frq arrays */
  rewind(fp);
  fscanf(fp, "%d", &i);
#ifndef NDEBUG
  assert(i = *nloci);
#endif
  while( fscanf(fp, "%d", &grp) != EOF)
  {
    i = findndx(grp, gid, *ngrps);
#ifndef NDEBUG
    assert(i >= 0);
    assert(i < *ngrps);
#endif
    for(j=0; j < *nloci; j++)  
    {
      if(fscanf(fp, "%d", &count) == EOF)
      {
	fprintf(stderr,"\nUnexpected EOF");
	exit(1);
      }
      if(count < 0)  /* ignore missing values */
	continue;
#ifndef NDEBUG
      assert(count <= *maxrepeat);
#endif
      (*frq)[i][j][count] += 1;
    }
  }
}







