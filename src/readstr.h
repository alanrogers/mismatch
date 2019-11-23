#ifndef READSTR_H
#define READSTR_H
/*** prototypes ***/
int findndx(int x, int *v, int len);
void readstr(
  FILE *fp,       /* pointer to input stream */
  int *gid,       /* gid[i] is identifier of i'th group */
  int *smp,       /* smp[i] = number of individuals in group i */
  int maxgroups,  /* dimension of arrays smp[] and gid[] */
  int ****frq,     /* (*frq)[i][j][k] = # w/ count k in locus j, group i */
  int *ngrps,      /* *ngrps is the number of groups */
  int *nloci,     /* *nloci is the number of loci */
  int *maxrepeat  /* *maxrepeat is the max repeat score */
 );
#endif /* READSTR_H */
