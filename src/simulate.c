#include <stdio.h>
#include <math.h>
#include "memstack.h"
#include "mismatch.h"
#include "bye.h"
#include "savetree.h"

/* externals */
extern int *mismatch, ***intermatch;            /* see iscoales.c */
extern  int nsubs, nwi, sampsize, n_sites;      /* ditto          */
int  total_iterations=0;    /* total number of iterations in all */
int  curr_iteration=0;      /* this iteration */


/* simulate a data set and estimate parameters */
int simulate(SIMULATION *s, POPHIST *history, double *theory_f, int msize,
	     MUTATION_MODEL mut_model, unsigned sflags)
{
  int mlen;
#ifdef SAVETREE
  FILE *fp;
  NODE *root;
#endif

  s->tree = iscoales(s->nsubs, s->subsiz, history, msize, 1.0);
#if 0
  /* mean tree depth is same as in version 2b */
  printf("\nsimulate: treedepth = %f", treedepth(s->tree));
#endif
  s->e.seg = getsegregating();

#ifdef SAVETREE
  if( treediff(s->tree, s->tree) )
  {
    printf("\ntreediff says tree != itself\n");
    exit(1);
  }else
  {
    printf("\ntreediff says tree == itself\n");
  }

  fp = mustopen("tree.bin", "w");
  writetree(s->tree, fp);
  fclose(fp);
  fp = mustopen("tree.bin", "r");
  root = readtree(fp);
  fclose(fp);

  if( treediff(s->tree, root) )
  {
    printf("\ntreediff says tree != copy\n");
    exit(1);
  }else
    printf("\ntreediff says tree == copy\n");
#endif

  if(s->tree==NULL)  
  {
    fprintf(stderr,"\nBad return from iscoales.");
    return(0);  /* error return */
  }
  mlen = getmatch(s->tree);
  estimate(s, msize, theory_f, mismatch, sflags);
  if(sflags & ROUGHNESS)
    s->e.roughness = getroughness(mismatch, msize);
  /* memory should be freed by calling routine */
#ifdef SAVETREE
  freetree(root);
#endif
  return(1);
}


/****************************************************************
Generate simulated data sets, and estimate parameters from each.
Resulting estimates are stored in array x.  The j'th col of x holds
estimates from the j'th iteration, ignoring iterations that don't
converge.  Do determine which variables go into which rows, consult
the indices defined in ndx.c.

An iteration is considered to have converged only if *all* requested
estimates are obtained.

The function returns the total number of iterations performed, including
those that did not converge.
****************************************************************/
int multisim(SIMULATION *sim, POPHIST *history, int iterations, double **x,
	     double *theory_f, int msize, MUTATION_MODEL mut_model,
	     unsigned sflags)
{
  int ntot, ngood, converged;
  extern int n_sites;
  extern int verbose;

  for(ngood=ntot=0; ngood < iterations; ntot++)
  {
    if(verbose)
    {
#if 0      
      putc('.', stderr);
#endif
      if(total_iterations > 0)
	fprintf(stderr,"\rprogress: %5.2f%%",
		100*(((double)curr_iteration++) / total_iterations));
      else
	fprintf(stderr,"\rit: %10d", curr_iteration++);
    }
    
    /* Perform simulation */
    if(simulate(sim, history, theory_f, msize, mut_model, sflags) == 0)
    {
      fflush(stdout);
#if 0      
      stackstatus(stderr);
#endif
      stackfree();  /* free memory used in simulation */
      sim->tree = NULL;
      continue;
    }
    /* copy estimates into column ngood of output array */
    converged = copy_estimates(x, ngood, &(sim->e), sflags);
#if 0    
    /* print estimates */
    printf("\nEstimates:\n");
    for(i=0; i<dim_s; i++)
      printf(" %10s", lbl_s[i]);
    putchar('\n');
    for(i=0; i<dim_s; i++)
      printf(" %10.6f", x[i][ngood]);
#endif
    if(converged)
      ngood++;
    fflush(stdout);
    stackfree();  /* free memory used in simulation */
    sim->tree = NULL;
  }
  return(ntot);
}
/*
 * Copy estimates from structure e into column "col" of array x.
 * Return 1 if all the estimates are good values, or 0 if any
 * are bad.
 */
int copy_estimates(double **x, int col, struct estimates *e, unsigned flags)
{
  int row, nrow=0, converged=1;

  if(flags & THETA0)
    x[nrow++][col] = e->theta0;
  if(flags & THETA1)
    x[nrow++][col] = e->theta1;
  if(flags & TAU)
    x[nrow++][col] = e->tau;
  if(flags & MSE)
    x[nrow++][col] = e->mse;
  if(flags & MAE)
    x[nrow++][col] = e->mae;
  if(flags & ROUGHNESS)
    x[nrow++][col] = e->roughness;
  if(flags & SEG)
    x[nrow++][col] = e->seg;

  for(row=0; row<nrow; row++)
  {
    if(x[row][col] == BADDBL)
    {
      converged = 0;
      break;
    }
  }
  return(converged);
}
