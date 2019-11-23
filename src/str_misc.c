/* MISCELLANEOUS FUNCTIONS FOR STRs */
#include <stdio.h>
#include "mytypes.h"
#include "iscoales.h"
#include "str_misc.h"

/* get STR frequency distribution from tree */
void str_getfrq(INTVEC *g, NODE *node)
{
  if(node == NULL)
    return;
  if(node->state.steps < g->lo)
    g->f[g->lo] += 1;
  else if(node->state.steps > g->hi)
    g->f[g->hi] += 1;
  else
    g->f[node->state.steps] += 1;
  str_getfrq(g, node->left);
  str_getfrq(g, node->right);
}

/* get STR mismatch distribution from frequency distribution */
void get_mismatch(int *mm, int *f1, int *f2, int maxrepeat)
{
  int i, j, diff;

  for(i=0; i<=maxrepeat; i++)
    mm[i] = 0;

  for(i=0; i<= maxrepeat; i++)
  {
    for(j=1; j<i; j++)
    {
      diff = i - j;
      mm[diff] += f1[i]*f2[j];
      mm[diff] += f1[j]*f2[i];
    }
    if(f1 != f2)
      mm[0] += f1[i]*f2[i];
    else
      mm[0] += f1[i]*(f1[i]-1);
  }
  /* mm[i] is now twice the correct value.  Fix it. */
  for(i=1; i<=maxrepeat; i++)
    mm[i] /= 2;
}

