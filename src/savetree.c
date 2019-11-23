#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mytypes.h"
#include "iscoales.h"
#include "bye.h"
#include "savetree.h"
extern int n_sites;
extern MUTATION_MODEL mut_model;

/* write a tree to file ofp */
int writetree(NODE *node, FILE *ofp)
{
  int rval = 0;
  
  if(node == NULL)
    return(0);

  /* write current node */
  if(fwrite(node, sizeof(NODE), 1, ofp) != 1) 
  {    
    fprintf(stderr,"\nwritetree[1]: bad rtn from fwrite\n");
    exit(1);
  }

  /* write list of descendants for current node */
  if(node->d->n > 0)
  {
    if(fwrite(node->d->ndx, sizeof(int), node->d->n, ofp) != node->d->n)
    {
      fprintf(stderr,"\nwritetree[2]: bad rtn from fwrite\n");
      exit(1);
    }
  }

  /* write genetic data for current node */
  if(node->state.seq != NULL)
  {
    assert(mut_model==finite_flat || mut_model==finite_gamma);
    assert(n_sites > 0);
    if(fwrite(node->state.seq, sizeof(SITE), n_sites, ofp) != n_sites)
    {
      fprintf(stderr,"\nwritetree[3]: bad rtn from fwrite\n");
      exit(1);
    }
  }

  rval = 1;
  rval += writetree(node->left, ofp);    /* write left branch */
  rval += writetree(node->right, ofp);   /* write right branch */
  
  return(rval);               /* return count of nodes written */
}

/* read a tree from file ifp */
NODE *readtree(FILE *ifp)
{
  NODE *node;

  node = (NODE *) malloc(sizeof(NODE));
  if(node == NULL)
  {
    fprintf(stderr,"\nreadtree: can't allocate node\n");
    exit(1);
  }

  /* read node from file */
  if( fread(node, sizeof(NODE), 1, ifp) != 1)
  {                              /* EOF */
    free(node);
    return(NULL);
  }

  /* read list of descendants for current node */
  if(node->d->n > 0)
  {
    node->d->ndx = (int *) mustalloc(node->d->n * sizeof(int));
    if(fread(node->d->ndx, sizeof(int), node->d->n, ifp) != node->d->n)
    {
      fprintf(stderr,"\nreadtree[2]: bad rtn from fread\n");
      exit(1);
    }
  }

  /* read genetic data for current node */
  if(node->state.seq != NULL)
  {
    assert(mut_model==finite_flat || mut_model==finite_gamma);
    assert(n_sites > 0);
    node->state.seq = (SITE *) malloc(n_sites * sizeof(SITE));
    if(node->state.seq == NULL)
    {
      fprintf(stderr,"\nreadtree: can't allocate state.seq\n");
      exit(1);
    }
    if(fread(node->state.seq, sizeof(SITE), n_sites, ifp) != n_sites)
    {
      fprintf(stderr,"\nreadtree[3]: bad rtn from fread\n");
      exit(1);
    }
  }

  /* connect left branch */
  if(node->left != NULL)
    node->left = readtree(ifp);

  /* connect right branch */
  if(node->right != NULL)
    node->right = readtree(ifp);
  return(node);
}

/* return 0 if trees are identical, 1 if they differ */
int treediff(NODE *n1, NODE *n2)
{
  int i;
  
  if(n1==NULL && n2 == NULL)
    return(0);

  if(n1 == NULL && n2 != NULL)
    return(1);

  if(n1 != NULL && n2 == NULL)
    return(1);

  if(n1->branch != n2->branch)
    return(1);

  if(n1->pop != n2->pop)
    return(1);

  if(n1->pop0 != n2->pop0)
    return(1);

  if(n1->state.steps != n2->state.steps)
    return(1);

  if(n1->state.seq!=NULL && n2->state.seq==NULL)
    return(1);

  if(n1->state.seq==NULL && n2->state.seq!=NULL)
    return(1);

  if(n1->state.seq!=NULL && n2->state.seq!=NULL)
  {
    assert(mut_model==finite_flat || mut_model==finite_gamma);
    assert(n_sites > 0);
    for(i=0; i<n_sites; i++)
    {
      if(n1->state.seq[i] != n2->state.seq[i])
	return(1);
    }
  }

  if(n1->mutations != n2->mutations)
    return(1);

  if(n1->d->n != n2->d->n)
    return(1);

  for(i=0; i<n1->d->n; i++)
  {
    if(n1->d->ndx[i] != n2->d->ndx[i])
      return(1);
  }

  if(treediff(n1->left, n2->left))
    return(1);

  if(treediff(n1->right, n2->right))
    return(1);

  return(0);
}

void freetree(NODE *node)
{
  if(node==NULL)
    return;

  freetree(node->left);
  freetree(node->right);
  if(node->state.seq != NULL)
    free(node->state.seq);
  free(node);
}


