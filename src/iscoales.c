#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "mytypes.h"
#include "bye.h"
#include "unirand.h"
#include "memstack.h"
#include "iscoales.h"  
#include "alloc2d.h"
#include "getcic.h"
#define INIT_SITEVAL 0      /* initial value of each site */
#define MUTATE(x) (!(x))    /* converts 0 to 1 and vice versa */
#define INBUFF 100
#define KEEPSIZE 100

/*** externals ***/
int *mismatch=NULL, ***intermatch=NULL;  /* set in alloc_match */
NODE *nodevec=NULL;                      /* ...... alloc_match */
int  **xtab = NULL;                      /* ...... alloc_match */
int msize=0;                             /* ...... alloc_match */
int nnodes=0, nextnode=0;                /* ...... iscoales    */
int n_sites=0;
int sampsize=0;     /* size of sample, all subdivisions included */
int nsubs=0;        /* number of subdivisions in sample */
int *subsize=NULL;  /* subsize[i] = size of i'th subdivision */
NODE ***smpl=NULL;   /* smpl[i][j] -> j'th sample from i'th subdivision*/

/** for mutations **/
double *cumprob=NULL;     /* cumprob[i] = prob that site is <= i */
char  *hit=NULL;        /* hit[i] = 1 (0) if site i has (has not) mutated*/
MUTATION_MODEL mut_model=infinite;  /* model of mutation */
int reset_mut=0;        /* should mutation rates be reset each time?*/
double gamma_shape;       /* gamma shape parameter */
double sim_mpd;           /* mean pairwise difference. Set by getmatch() */

/*--Allocate & free 2 arrays--*/
void  **stackalloc2d(unsigned dim1, unsigned dim2, unsigned size)
{
  int     i;
  unsigned nelem;
  char   *p, **pp;

  nelem = dim1 * dim2;
  p = (void *) stackalloc( nelem * size );
  if (p == NULL)
    return (NULL);
  pp = (char **) stackalloc( dim1 * sizeof(char *) );
  if (pp == NULL)
    return (NULL);
  for (i = 0; i < dim1; i++)
    pp[i] = p + i * dim2 * size;

  return ((void **) pp);
}
/* Join two nodes into a single node */
NODE   *newnode(NODE * n1, NODE * n2)
{
  NODE   *new;

#ifndef NDEBUG
  if(n1!=NULL && n2!=NULL && n1->pop != n2->pop)
  {
    fprintf(stderr,
	    "\nError in newnode: %d = n1->pop != n2->pop=%d\n",
	    n1->pop, n2->pop);
    exit(1);
  }
#endif  
  new = (NODE *) nodevec + nextnode;
  ++nextnode;
  new->branch = 0.0;
  new->left = n1;
  new->right = n2;
  new->d = NULL;
  new->ancestor = NULL;
  new->mutations = 0;
  switch(mut_model)
  {
      case finite_flat:
      case finite_gamma:
	  new->state.seq = NULL;
	  break;
      case stepwise:
	  new->state.steps = 0;
	  break;
      default:
	  break;
  }
  if(n1 != NULL)
    n1->ancestor = new;
  if(n2 != NULL)
    n2->ancestor = new;
  if(n1 != NULL)
    new->pop = n1->pop;
  else if(n2 != NULL)
    new->pop = n2->pop;
  else
    new->pop = -1;
  if(n1 != NULL)
    new->d = n1->d;
  if(n2 != NULL)
    new->d = mergelist(new->d, n2->d);
  return (new);
}
#ifndef NDEBUG
void check_s(char *msg, int *s, int K, NODE **node, int S, int SS)
{
  int sum, sum2, count[MAX], i, bad=0;

  for(i=0; i<K; i++)
    count[i] = 0;
  for(i=0; i<S; i++)
  {
    if((node[i]->pop) < 0 || (node[i]->pop) >=K)
    {
      fprintf(stderr,"\ncheck_s:node[%d]->pop=%d. Legal=0..%d",
	      i, node[i]->pop, K-1);
      bad = 1;
    }
    count[node[i]->pop] += 1;
  }
  for(sum=sum2=i=0; i<K; i++)
  {
    sum += count[i];
    sum2 += count[i]*count[i];
  }
  if(sum != S)
  {
    fprintf(stderr,"\ncheck_s: sum=%d != %d=S", sum, S);
    bad = 1;
  }
  if(sum2 != SS)
  {
    fprintf(stderr,"\ncheck_s: sum2=%d != %d=SS", sum2, SS);
    bad = 1;
  }
  for(i=0; i<K; i++)
  {
    if(count[i] != s[i])
    {
      fprintf(stderr,"\ncheck_s: count[%d]=%d!=%d=s[%d]",
	      i, count[i], s[i], i);
      bad = 1;
    }
  }
  if(bad)
  {
    fprintf(stderr,"\ncheck_s: %s, S=%d\ns=", msg, S);
    for(i=0; i<K; i++)
      fprintf(stderr," %d", s[i]);
    exit(1);
  }
}
#endif

double   getr(double S, double SS, double mn, double theta)
{
  return( ( S*mn + 0.5*(SS - S) )/theta );
}

/****************************************************************
Build a random tree using coalescent w/ geographic structure and
island model migration.

Let s[i] = # of nodes in pop i
    S = sum(s[i]) = total nodes
    SS = sum(s[i]*s[i]), the sum of squares
    n = size of a subpopulation
    K = # of subpopulations
    m = island model migration rate
    a = probability density that an event will happen in (t, t+dt)
        (i.e. the hazard)

   K-1
   ---
   \
a = >   (  s[i]*m + (s[i]*(s[i]-1)/2) / n )
   /__
   i=0

  = S*m + ( SS - S)/(2*n)

I will only need a in calculating a*t, where t is time in generations.
This is equal to

      2*u*t
a*t = ----- * [ S*(mn) + (SS - S)/2 ]
      2*u*n

Thus, if I measure time by tau = 2*u*t, then the parameters are:
theta=2*u*n and mn = m*n.  In mutational time, the rate a becomes

r = [ S*(mn) + (SS - S)/2 ]/theta

The time of the next event is an exponential r.v. with rate r.  The
cdf is 1 - exp(-r*x).  Since any cdf is uniformly distributed on
(0,1), one can generate the exponential variate as follows.  First
generate a uniform variate u, and then solve u = 1 - exp(-r*x) to get 
x = -ln(1-u)/r.

Given that an event has occurred, it is a migration with probability 

 pmig = S*(mn)/[ S*(mn) + (SS - S)/2]

When a coalescent event occurs in group i, the following adjustments
are made:
1. 2 nodes are joined to form a single new node
2. S -= 1;
3. SS -= 2*s[i] - 1;
4. s[i] -= 1;
Here, step 3 is equivalent to

  SS' = SS - s[i]^2 + (s[i] - 1)^2

When a migration from i to j occurs, the following adjustments are
made:
1. SS += 2*(s[j] - s[i] + 1);
2. s[i] -= 1;
3. s[j] += 1;
Here, step 1 is equivalent to

  SS' = SS - s[i]^2 + (s[i]-1)^2 - s[j]^2 + (s[j]+1)^2

ON ENTRY:
  sampsize = # of nodes per subpopulation
  mut_rate: should equal 1.0 to scale time by mutation rate.
****************************************************************/
NODE   *iscoales(int init_nsubs, int *init_subsize, POPHIST *history,
		 int init_msize, double mut_rate)
{
  POPHIST *ph;
  NODE   **node;
  int     s[MAX];  /*s[i]=#of nodes in pop i */
  int     S, SS;
  int     i, j, k, ii, jj, ipop, jpop, oldK, kk, ndx;
  int     v[MAX];
  double   t_g;     /*time since last migration or coalescence */
  double   t_p;     /*time since last change in population parameters */
  double   c[MAX];  /* cumulative distribution function */
  double   pmig;
  double   x, r;

  clear_mutation();

  /****************************************************************
  Set permanent arrays, which do *not* get cleared each time.  If
  nsubs and msize have not changed, I assume w/o checking that
  subdivision sizes haven't changed either.
  ****************************************************************/
  if(msize != init_msize || nsubs!=init_nsubs)
  {
    clear_externals();
    msize = init_msize;
    nsubs = init_nsubs;
    subsize = (int *) mustalloc(nsubs * sizeof(int));
    sampsize=0;
    for(i=0; i<nsubs; i++)
      sampsize += subsize[i] = init_subsize[i];
    nnodes = 2*sampsize - 1;      /* total nodes in tree */
    mismatch = (int *) mustalloc(msize * sizeof(int));
    /* nodevec: pointers that WILL get reshuffled */
    nodevec = (NODE *) mustalloc(nnodes * sizeof(NODE));
    xtab = (int **) alloclt(sampsize, sizeof(int));
    if(xtab==NULL)
      error("iscoales: memory");
    intermatch = (int ***) alloc3d(nsubs, nsubs, msize, sizeof(int));
    if(intermatch == NULL)
	error("iscoales: memory");

    /* smpl: pointers that will not be reshuffled by coalescent events */
    smpl = (NODE ***) mustalloc(nsubs * sizeof(NODE **));
    for(i=0; i<nsubs; i++)
      smpl[i] = (NODE **) mustalloc(subsize[i] * sizeof(NODE *));
  }

  /****************************************************************
  The number of subdivisions in the sample may not match the number
  in the population.  If not, I add another epoch to the population
  history with length 0.  The values of mn and theta don't matter
  since this epoch has 0 length anyway.
  ****************************************************************/
  if(history->K != nsubs)
  {
    ph = newhistory(1.0, 1.0, 0.0, nsubs); /*new epoch*/
    ph->next = history;                    /*attach old history*/
    history = ph;
  }

  /****************************************************************
  The newnode function doesn't allocate memory; it just grabs the
  next available node from the nodevec vector.  Setting nextnode=0
  tells newnode to start at the beginning again, effectively freeing
  all the nodes that were used last time.
  ****************************************************************/
  nextnode = 0;
  node = (NODE **) stackmustalloc(sampsize * sizeof(NODE*));
  for(k=j=0; j<nsubs; j++)
  {
    for (i = 0; i < subsize[j]; i++)
    {
      smpl[j][i] = node[k] = newnode(NULL, NULL);
      node[k]->d = newlist(1, &k);
      node[k]->pop = j;
      node[k]->pop0 = j;
      k += 1;
    }
    s[j] = subsize[j];
  }
  /****************************************************************
  At this point we can access the nodes in the sample either through
  smpl[j][i] or through node[k].  After the building the tree, the
  entries of node[] will will be scrambled but those of smpl[][] will
  still point to the sample.
  ****************************************************************/

  S  = sampsize;                  /* number of nodes */
  for(SS=0.0, i=0; i<nsubs; i++)
    SS += subsize[i]*subsize[i];  /*sum of squared sample sizes */
  if (S > MAX)
  {
    fprintf(stderr, "\nError: trees can only have %d tips\n", MAX);
    exit(1);
  }
#ifndef NDEBUG    
  assert(S == sampsize);
  check_s("top of iscoales", s, history->K, node, S, SS);
#endif
  if(history->tau==0.0)
    history->theta=1.0;

  /* initial value of rate r of exponential random variable */
  r = getr((double) S, (double) SS, history->mn, history->theta);

  /* Iterate until sample size (S) equals 1 */
  t_g = t_p = 0.0;
  while(S > 1)
  {
    x = -log(1.0-uni());  /* An exponential random variate. */
    /* translate x into units of 1/(2u) generations */
    while(history->next != NULL  && (r == 0 || t_p + x/r > history->tau))
    {
      /****************************************************************
      1st line below subtracts off portion of x that was "used up" by the
      epoch we are about to leave.
      ****************************************************************/
      x -= r*(history->tau - t_p);
      /****************************************************************
      The rest of this loop resets parameters to reflect the new epoch
      and recalculates r.
      ****************************************************************/
      t_g += history->tau - t_p;
      t_p = 0.0;
      oldK = history->K;
      history = history->next;
      if(history->tau==0.0)
	history->theta=1.0;
      if(history->K < oldK)  /*  reduce # of subpopulations */
      {
	/* get vector to map oldK old pops into history->K new ones */
	get_collapse_vector(oldK, history->K, v);
	for(i=0; i<oldK; i++)
	  s[i] = 0;
	for(i=0; i<S; i++)
	{
	  node[i]->pop = v[node[i]->pop];  /* reset population values */
	  s[node[i]->pop]++;               /* recalc pop sample sizes */
	}
	SS = 0;                            /* reset sum of squares, SS*/
	for(i=0; i<history->K; i++)
	  SS += s[i]*s[i];
      }else if(history->K > oldK)  /* increase # of supopopulations */
      {
	/***************************************************
	Initially, the new groups will be empty: s[i]=0.
	These groups gain members only by migration.   
	This assumes that, in forwards time, the reduction in group
	numbers was caused by group extinction rather than
	coalescence. 
        ****************************************************/
	for(i=oldK; i<history->K; i++)
	  s[i] = 0;
      }
      /* reset r for new history parameters */
      r = getr((double) S, (double) SS, history->mn, history->theta);
    }
    if(r == 0)
      error("iscoales: NO COALESCENCE\n");
    t_p += x/r;  /* time since last change in history parameters */
    t_g += x/r;  /* time since last coalescence */

#ifndef NDEBUG    
    check_s("after setting t_p & t_g", s, history->K, node, S, SS);
#endif

    /* Classify current event */
    pmig = (double) S*(history->mn);
#ifndef NDEBUG    
    assert(pmig + 0.5*(SS-S) > 0.0);
#endif
    pmig = pmig / (pmig + 0.5*(SS-S));  /* Pr[event is a migration] */
    if(history->mn > 0 && uni() < pmig)
    {
      /* event is a migration: choose migrant */
      i = randint(S);
      ipop = node[i]->pop;
#ifndef NDEBUG
      assert(s[ipop] > 0);
#endif
      /* choose group of origin */
      jpop = randint((history->K)-1);
      if(jpop >= ipop)
	jpop++;
#ifndef NDEBUG
      assert(ipop<history->K);
      assert(jpop<history->K);
#endif
      /* move node i */
      node[i]->pop = jpop;
      /* reset parameters */
      SS += 2*(s[jpop] - s[ipop] + 1);
      s[ipop] -= 1;
      s[jpop] += 1;
#if 0
    check_s("after migration", s, history->K, node, S, SS);
#endif
      continue;
    }else /* Event is a coalescence. */
    {
      /* record time */
      for (i = 0; i < S; i++)
	node[i]->branch += t_g;
      t_g = 0.0;
#ifndef NDEBUG
      assert(SS > S);
#endif
      /* Choose group */
      x = SS - S;
      c[0] = s[0]*(s[0]-1)/x;
      for(i=1; i<history->K; i++)
	c[i] = c[i-1] + s[i]*(s[i]-1)/x;
      c[history->K - 1] = 1.0;
      x = uni();
      for(ipop=0; ipop<history->K && c[ipop]<=x; ipop++)
	;
      /* Choose two nodes w/i group ipop */
#ifndef NDEBUG
      assert(s[ipop] >= 2);
#endif
      i = randint(s[ipop]);
      j = randint(s[ipop] - 1);
      if(j >= i)
	j += 1;
      /* find nodes i & j w/i group ipop */
      ndx = ii = jj = -1;
      for(kk=0; kk<S; kk++)
      {
	if(node[kk]->pop == ipop)
	{
	  ndx++;
	  if(ndx == i)
	    ii = kk;
	  else if(ndx == j)
	    jj = kk;
	  if(ii>=0 && jj>=0)
	    break;
	}
      }
#ifndef NDEBUG
      if(ii<0 || jj<0)
      {
	fprintf(stderr,"\nbad ii or jj: ii=%d jj=%d ipop=%d", ii, jj, ipop);
	fprintf(stderr,"  s[%d]=%d", ipop, s[ipop]);
	fprintf(stderr,"\n   node[*]->pop=");
	for(i=0; i<S; i++)
	{
	  fprintf(stderr," %d", node[i]->pop);
	  if((i+1)%20==0)
	    putc('\n', stderr);
	}
	exit(1);
      }
      if(node[ii]->pop != node[jj]->pop)
      {
	fprintf(stderr,
		"\niscoales: node[%d]->pop=%d != %d=node[%d]->pop",
		ii,node[ii]->pop, node[jj]->pop, jj);
	fprintf(stderr," S=%d\n", S);
	exit(1);
      }
#endif
      /* join the two nodes */
      node[ii] = newnode(node[ii], node[jj]);
      /* turn the S-vector into an (S-1)-vector */
      if(ii == S-1)
	node[jj] = node[ii];
      else if(jj != S-1)
	node[jj] = node[S-1];
      /* reset parameters */
      S -= 1;
      SS -= 2*s[ipop] - 1;
      s[ipop] -= 1;
    }
    /* reset r because SS and maybe S have changed */
    r = getr((double) S, (double) SS, history->mn, history->theta);
  }
  mutate(node[0], NULL, mut_rate);  /* add mutations to tree */
  return(node[0]);
}
/****************************************************************
Read file fp, which should contain one line of data for each time
period, with the most recent time period first.  Each line should
contain four items of data: theta, mn, tau, and K.  Here theta=4*u*n*K
where n is the size of a single population.

The function returns a pointer to the base of a linked list of type
POPHIST. 
****************************************************************/
POPHIST *gethistory(FILE *fp)
{
  double mn, theta, tau;
  int i, K, got_infinity = 0;
  POPHIST *root=NULL, *ph, *last=NULL;
  char buff[INBUFF];

  while( getdoubleic(&theta, buff, INBUFF, fp)!=EOF)
  {
    if(got_infinity)
      error("pophist.ini continues after an infinite epoch");
    if( getdoubleic(&mn, buff, INBUFF, fp)==EOF )
      error("Unexpected EOF reading mn in pophist.ini");
    if( getwordic(buff, INBUFF, fp)==NULL )
      error("Unexpected EOF reading tau in pophist.ini");
    /* check for word "inf" */
    for(i=0; i<INBUFF && buff[i] != '\0'; i++)
      buff[i] = tolower(buff[i]);
    if( strncmp(buff, "inf", INBUFF-1)==0 )
    {                   /* doesn't matter what number we put here */
      tau = 99.99;      /* because initial value of tau is ignored*/
      got_infinity = 1;
    }else if(!isdigit(*buff))
    {
      fprintf(stderr,"\nIllegal value of tau: \"%s\"\n", buff);
      exit(1);
    }else
      sscanf(buff,"%lf", &tau);
    if( getintic(&K, buff, INBUFF, fp)==EOF )
      error("Unexpected EOF reading K in pophist.ini");
    
    if(K < 1)
    {
      fprintf(stderr,"\nerror in population history data: K=%d\n", K);
      exit(1);
    }
    ph = newhistory(theta/K, mn, tau, K);
    if(last==NULL)
      root = ph;        /* set root on first pass through loop */
    else
      last->next = ph;  /* append to list on later passes */
    last = ph;          /* always set last to end of list */
  }
  if( !got_infinity )
  {
    fprintf(stderr,
     "\nWarning: tau was not infinite in initial epoch of pophist.int.");
    fprintf(stderr,"\n         I'm treating it as infinite anyway.");
  }
  return(root);
}

/*** create a new structure of type POPHIST ****/
POPHIST *newhistory(double theta, double mn, double tau, int K)
{
  POPHIST *ph;
  
  ph = (POPHIST *) mustalloc( sizeof(POPHIST));
  ph->next = NULL;
  ph->mn = mn;
  ph->theta = theta;
  ph->tau = tau;
  ph->K = K;
  if(K==1)
    ph->mn = 0.0;   /* can't have migration w/ only 1 group */
  return(ph);
}
/*** free a list of type POPHIST ***/
void freehistory(POPHIST *ph)
{
  if(ph == NULL)
    return;
  freehistory(ph->next);
  free(ph);
}

/*** duplicate a list of type POPHIST ****/
POPHIST *duphist(POPHIST *new, POPHIST *old)
{
  if(old == NULL)
  {
    freehistory(new);
    return(NULL);
  }
  if(new==NULL)
  {
    new = newhistory(old->theta, old->mn, old->tau, old->K);
    new->next = NULL;
  }else
  {
    new->theta = old->theta;
    new->mn = old->mn;
    new->tau = old->tau;
    new->K = old->K;
  }
  new->next = duphist(new->next, old->next);
  return(new);
}
void writehistory(FILE *fp, POPHIST *history, int comment)
{
  fprintf(fp,"\n%c%10s %10s %10s %10s", comment, "theta", "mn", "tau", "K");
  if(history==NULL)
    return;
  while(history->next != NULL)
  {
    fprintf(fp, "\n%c%10.4f %10.4f %10.4f %10d",
	    comment,
	    (history->theta) * (history->K),
	    history->mn,
	    history->tau,
	    history->K);
    history = history->next;
  }
  fprintf(fp, "\n%c%10.4f %10.4f %10s %10d",
	  comment,
	  (history->theta) * (history->K),
	  history->mn,
	  "Inf",
	  history->K);
  putc('\n', fp);
}
/********integer comparison function used by qsort()************/
int icompar(int *x, int *y)  
{
  if(*x > *y)
    return(1);
  if(*x < *y)
    return(-1);
  return(0);
}
/****************************************************************
The collapse_vector, v, is used to reduce the number of groups from K1
to K2.  Each individual in old pop i is assigned to new pop v[i].  To
generate v, I calculate a random partition of the old groups.
****************************************************************/
int get_collapse_vector(int K1, int K2, int *v)
{
  int i, j, r[MAX], p[MAX];

  (void) randperm(r, K1);    /* randomly reorders groups */
  (void) randperm(p, K1-1);  /* vector of random partition points */
  if(K2 > 2)                 /* sort 1st K2-1 partition points */
    qsort(p, (unsigned) (K2-1), sizeof(int),
     (int (*)(const void*,const void*)) icompar);
  for(i=j=0; i<K2-1; i++)   /* j indexes partitions */
  {
    while(j <= p[i])        /* we are in the i'th partition */
      v[r[j++]] = i;        /* assign i to all members of this partition */
  }                         /* r[] randomizes partition membership */
  while(j<K1)               /* the last partition goes to the end */
    v[r[j++]] = i;
  return(0);
}

void    prtree(NODE * node, int indent)
{
  int     i;

  if (node == NULL)
    return;
  putchar('\n');
  for (i = 0; i < indent; i++)
    putchar('-');
  printf("lngth=%f pop=%d mut=%d",node->branch, node->pop, node->mutations);
  switch(mut_model)
  {
      case stepwise:
	  printf(" steps=%d", node->state.steps);
	  break;
      case finite_gamma:
      case finite_flat:
	  putchar(':');
	  for(i=0; i<n_sites; i++)
	      printf("%c", node->state.seq[i]);
	  break;
      default:
	  break;
  }
  prtree(node->left, indent + 2);
  prtree(node->right, indent + 2);
}

/***** measure depth of tree *****/
double treedepth(NODE *node)
{
  if(node == NULL)
    return(0.0);
  return(node->branch + treedepth(node->left));
}

/**** measure length of tree = sum of all branch lengths ****/
double treelength(NODE *node)
{
  if(node == NULL)
    return(0.0);
  return(node->branch + treelength(node->left) + treelength(node->right));
}

/**** count mutations in tree ****/
int treemutations(NODE *node)
{
  if(node == NULL)
    return(0);
  return(node->mutations + treemutations(node->left)
	 + treemutations(node->right));
}
/**** Get mismatch and intermatch distributions ******************/
int getmatch(NODE *tree)
{
  register int dif, sum;
  int i, j, ipop, jpop;
  int size;

#ifndef NDEBUG
  assert(sampsize > 0);
  assert(msize > 0);    /* make sure externals have been initialized */
  assert(mut_model==infinite || mut_model==finite_flat
	 || mut_model==finite_gamma || mut_model == stepwise);
#endif  

  memset((void *) mismatch, 0, sizeof(int)*msize);   /* initialize */
  memset((void *) (intermatch[0][0]), 0, sizeof(int)*nsubs*nsubs*msize);

  sim_mpd = 0.0;
  switch(mut_model)
  {
  case finite_flat:
  case finite_gamma:
    for(i=1; i < sampsize; i++)
    {
#ifndef NDEBUG      
      assert(nodevec[i].pop0 >= 0);
#endif
      for(j=0; j<i; j++)
      {
	sim_mpd += dif = ndiffs(nodevec+i, nodevec+j);
	if(dif >= msize)
	  dif = msize-1;             /* dump extreme vals into last entry */
	++mismatch[dif];             /* accumulate mismatch counts */
#ifndef NDEBUG
	assert(nodevec[j].pop0 >= 0);
#endif
	if(nodevec[i].pop0 > nodevec[j].pop0)
	{
	  ipop = nodevec[i].pop0;
	  jpop = nodevec[j].pop0;
	}else
	{
	  jpop = nodevec[i].pop0;
	  ipop = nodevec[j].pop0;
	}
	++intermatch[ipop][jpop][dif];
      }
    }
    break;
  case stepwise:
  case infinite:
    /****************************************************************
    fill xtab w/ zeroes: Note that we have
    sampsize*(sampsize*(sampsize+1))/2 entries rather than
    sampsize*(sampsize*(sampsize-1))/2.  This is because alloclt() allocates
    a full lower triangular matrix (including the diagonal) even though
    the diagonal entries are not used.  If I were a better person I would
    eliminate this inefficiency.
    ****************************************************************/  
    memset((void *) xtab[0], 0, sizeof(int)*(sampsize*(sampsize+1))/2);
    /* put counts in xtab matrix: */
    crosstab(tree);
    /****************************************************************
    use xtab matrix to calculate mismatch distribution and mean
    pairwise difference
    ****************************************************************/
    for (i = 1; i < sampsize; i++)
    {
      for (j = 0; j < i; j++)
      {
	fflush(stdout);
	sim_mpd += dif = xtab[i][j];
	if(dif >= msize)
	  dif = msize-1;
	++mismatch[dif];
	if(nodevec[i].pop0 > nodevec[j].pop0)
	{
	  ipop = nodevec[i].pop0;
	  jpop = nodevec[j].pop0;
	}else
	{
	  jpop = nodevec[i].pop0;
	  ipop = nodevec[j].pop0;
	}
	++intermatch[ipop][jpop][dif];
      }
    }
    break;
#ifndef NDEBUG    
  default:
    fprintf(stderr,"\nIllegal switch value in getmatch\n");
    exit(1);
#endif
  }
  sim_mpd /= (sampsize*(sampsize-1))/2;
  /* discard trailing zeroes */
  for(size=msize; size>0 && mismatch[size-1]==0; --size)
    ;
    
#if 1  
  /* check mismatch distribution */
  for(sum=i=0; i<size; i++)
    sum += mismatch[i];
  if(sum!= (sampsize*(sampsize-1))/2)
    printf("\nmismatch:sum(m)=%d, should = %d",
	   sum, (sampsize*(sampsize-1))/2);
  for(ipop=0; ipop<nsubs; ipop++)
    for(jpop=0; jpop<=ipop; jpop++)
    {
      for(sum = i = 0; i<size; i++)
	sum += intermatch[ipop][jpop][i];
      if(ipop==jpop)
      {
	if(sum!=(subsize[ipop]*(subsize[ipop]-1)/2))
	{
	  fflush(stdout);
	  fprintf(stderr,"\nw/i group mismatch sum=%d; should be %d\n",
		  sum, (subsize[ipop]*(subsize[ipop]-1)/2));
	  exit(1);
	}
      }else
      {
	if(sum!=subsize[ipop]*subsize[jpop])
	{
	  fflush(stdout);
	  fprintf(stderr,"\nw/i group mismatch sum=%d; should be %d\n",
		  sum, (subsize[ipop]*subsize[jpop]));
	  exit(1);
	}
      }
    }
#endif
  return(size);
}

/* set external defs back to original state */
void clear_externals(void)
{
  int i;

  if(mismatch != NULL)
  {
    free(mismatch);
    mismatch = NULL;
  }
  if(intermatch != NULL)
  {
    free3d((void ***) intermatch);
    intermatch = NULL;
  }
  if(xtab != NULL)
  {
    free2d((void **) xtab);
    xtab = NULL;
  }
  if(nodevec!= NULL)
  {
    free(nodevec);
    nodevec = NULL;
  }
  if(smpl != NULL)
  {
    for(i=0; i<nsubs; i++)
      free(smpl[i]);
    free(smpl);
  }
  mismatch = (int *) NULL;
  intermatch = (int ***) NULL;
  xtab = (int **) NULL;
  nodevec = (NODE *) NULL;
  smpl = (NODE ***) NULL;
  sampsize = msize = nsubs = nnodes = nextnode = 0;
}

/****************************************************************
 Create a matrix whose ij'th entry is the difference between individuals
 i and j.
 ****************************************************************/
void    crosstab(NODE * node)
{
  int     nd, *descendant, nc, complement[MAX];
  int increment=0;
  register int i, j;

  if (node == NULL)
    return;

  switch(mut_model)
  {
  case stepwise:
    increment = node->state.steps;
    break;
  case infinite:
    increment = node->mutations;
    break;
  default:
    error("Wrong mutation model in crosstab");
  }
      
  if(increment > 0)
  {
    nd = node->d->n;
    descendant = node->d->ndx;
    /* Get complement of list of descendants */
    for (nc = i = j = 0; i < nd; i++)
    {
      while (j < descendant[i])
	complement[nc++] = j++;
      j++;
    }
    while (j < sampsize)
      complement[nc++] = j++;
    assert(nc + nd == sampsize);

    for (i = 0; i < nd; i++)
      for (j = 0; j < nc; j++)
      {
	if (descendant[i] > complement[j])
	  xtab[descendant[i]][complement[j]] += increment;
	else
	  xtab[complement[j]][descendant[i]] += increment;
      }
  }
  crosstab(node->left);
  crosstab(node->right);
}

LIST   *newlist(int n, int *ndx)
{
  int     i;
  LIST   *l;

  l = (LIST *) stackmustalloc( sizeof(LIST));
  l->n = n;
  l->ndx = (int *) stackmustalloc(n * sizeof(int));
  for (i = 0; i < n; i++)
    l->ndx[i] = ndx[i];
  return (l);
}

/* print list */
void prlist(LIST *l)
{
  int i;

  for(i=0; i<l->n; i++)
    printf(" %d", l->ndx[i]);
}

/* merge two lists */
LIST   *mergelist(LIST * l1, LIST * l2)
{
  LIST   *new;
  int     i1, i2, inew;

  if (l2 == NULL)
    return (l1);
  if (l1 == NULL)
    return (l2);

  new = (LIST *) stackmustalloc( sizeof(LIST) );
  new->n = l1->n + l2->n;
  new->ndx = (int *) stackmustalloc( new->n * sizeof(int) );
  i1 = i2 = inew = 0;
  while ((i1 < l1->n) && (i2 < l2->n))
  {
    if (l1->ndx[i1] < l2->ndx[i2])
      new->ndx[inew++] = l1->ndx[i1++];
    else
      new->ndx[inew++] = l2->ndx[i2++];
  }
  while (i1 < l1->n)
  {
    new->ndx[inew++] = l1->ndx[i1++];
  }
  while (i2 < l2->n)
  {
    new->ndx[inew++] = l2->ndx[i2++];
  }
#ifndef NDEBUG
  assert(i1 == l1->n && i2 == l2->n && inew == new->n);
#endif

  return (new);
}

/****************************************************************
init_mutation: call this to initialize mutation model
Parameters
 init_mut_model: specifies mutation model
 init_n_sites  : the number of sites.
 init_shape    : shape parameter of gamma distribution
 init_reset    : if nonzero, gamma rates are reset by each call to iscoales 
init_n_sites and init_shape are both ignored when mut_model==infinite.
****************************************************************/
double init_mutation(MUTATION_MODEL init_mut_model, int init_n_sites,
		   double init_shape, int init_reset)
{
  /** if mutation vectors are already set then unset them **/
  if(n_sites>0 && n_sites != init_n_sites)
  {
    if(cumprob != NULL)
    {
      free(cumprob);
      cumprob=NULL;
    }
    if(hit != NULL)
    {
      free(hit);
      hit = NULL;
    }
  }

  mut_model = init_mut_model;
  switch(mut_model)
  {
  case finite_gamma:
    reset_mut = init_reset;
    gamma_shape = init_shape;
    cumprob = (double *) mustalloc(init_n_sites * sizeof(double));
    n_sites = init_n_sites;
    hit = (char *) mustalloc(n_sites * sizeof(char));
    /****************************************************************
    If reset_mut is nonzero, gamma rates will be set at top of iscoales.
    Otherwise, set them here.
    ****************************************************************/
    if(!reset_mut)
      set_gamma_rates();
    break;
  case finite_flat:
    n_sites = init_n_sites;
    hit = (char *) mustalloc(n_sites * sizeof(char));
    break;
  case stepwise:
  case infinite:
    n_sites = 0;
    break;
  default:
    fflush(stdout);
    fprintf(stderr,"\ninit_mutation: illegal value of mut_model");
    exit(1);
  }
  return(0.0);  /* meaningless return value */
}

/****************************************************************
SET GAMMA MUTATION RATES:
The mutational time scale requires that the sum of mutation rates
across sites equal unity so that there will on average be one mutation
per unit of mutational time.  To accomplish this goal, we first define
$y\equiv x/\beta$, where $x$ is gamma-distributed with scale parameter
$\beta$ and shape parameter $\alpha$.  The variable $y$ is also
gamma-distributed, with density 
\begin{displaymath} 
      y^{\alpha-1} e^{-y} / \Gamma(\alpha) 
\end{displaymath} 
We first use this density to generate a vector
$(y_1,y_2,\ldots{},y_K)$.  These variates are each equal to $x/\beta$,
where $x$ is the gamma-distributed variate we seek and $\beta$ is
arbitrary.  To satisfy the constraint imposed by the mutational time
scale, we set $\beta = 1/\sum_i y_i$.  In other words, we set $x_i =
y_i / \sum_i y_i$.  This provides a vector of gamma-distributed
mutation rates that sums to unity, as required.
****************************************************************/
void set_gamma_rates(void)
{
  int i;
  double b;
  
  b = 0.0;
  for(i=0; i<n_sites; i++)
    b += cumprob[i] = gamma_dev((double) gamma_shape);
  cumprob[0] /= b;
  for(i=1; i<n_sites; i++)
    cumprob[i] = cumprob[i-1] + cumprob[i]/b;
}

/**** call before building each tree ****/
void clear_mutation(void)
{
  sim_mpd = -1.0;  /* indicates that sim_mpd has not been calculated */
  switch(mut_model)
  {
  case finite_gamma:
    if(cumprob==NULL)
    {
      fprintf(stderr,"\nerror: call init_mutation before clear_mutation\n");
      exit(1);
    }
    if(reset_mut)  /* set gamma mutation rates on each call */
      set_gamma_rates();
    /** fall through to finite_flat **/
  case finite_flat:
    if(hit==NULL)
    {
      fprintf(stderr,"\nerror: call init_mutation before clear_mutation\n");
      exit(1);
    }
    memset(hit, 0, (unsigned) n_sites);  /* set to zero */
    break;
  case stepwise:
  case infinite:
    n_sites = 0;
    break;
#ifndef NDEBUG
  default:
    fflush(stdout);
    fprintf(stderr,"\nclear_mutation: illegal value of mut_model");
    exit(1);
#endif
  }
}

/* put mutations onto tree */
void mutate(NODE *node, STATE *inherited, double mut_rate)
{
  int site, nmut;

  if(node==NULL)
    return;
  /* # of mutations */
  node->mutations = nmut = poidev(0.5*node->branch*mut_rate);

  switch(mut_model)
  {
  case stepwise:
    if(inherited == NULL)
      node->state.steps = getsteps(node->mutations);
    else
      node->state.steps = inherited->steps + getsteps(node->mutations);
    break;
  case infinite:
    n_sites += nmut;
    break;
  case finite_gamma:
  case finite_flat:
    if(node->mutations == 0 && inherited == NULL)
      node->state.seq = NULL;
    else
      node->state.seq = copyseq(inherited->seq);
    while(nmut-- > 0)
    {
      /* choose a site to mutate */
      site = getsite();
      hit[site] = 1;
      node->state.seq[site] = MUTATE(node->state.seq[site]);
    }
    break;
  default:
    error("mutate: Illegal mutation model");
  }
  mutate(node->left, &(node->state), mut_rate);
  mutate(node->right, &(node->state), mut_rate);
}

/* choose a site to mutate under one of the finite sites models */
int getsite(void)
{
  double r;
  int site, lo, mid, hi;

  switch(mut_model)
  {
  case finite_flat:
    site = randint(n_sites);
    break;
  case finite_gamma:
    /* choose a site to mutate by binary search */
    r = uni();  /* uniform random variable */
    lo = 0;
    hi = n_sites - 1;
    while(lo < hi)
    {
      mid = lo + (hi-lo)/2;
      if(mid == lo)
      {
	if(cumprob[lo]<r)
	  lo = hi;
	break;
      }
      if(cumprob[mid] <= r)
	lo = mid;
      else
	hi = mid;
    }
    return(lo);
    break;
  default:
    fflush(stdout);
    fprintf(stderr,"\nIllegal mutation model in getsite()\n");
    exit(1);
  }
  return(site);
}

/* copy DNA (or restriction site) sequence */
SITE *copyseq(SITE *s1)
{
  int i;
  SITE *s;

  s = (SITE *) stackalloc( n_sites*sizeof(SITE) );
  if(s==NULL)
  {
    fflush(stdout);
    fprintf(stderr,"\ncopyseq: stack full\n");
    exit(1);
  }
  if(s1 == NULL)
    for(i=0; i<n_sites; i++)
      s[i] = INIT_SITEVAL;
  else
    for(i=0; i<n_sites; i++)
      s[i] = s1[i];
  return(s);
}

/* count diffs between two nodes */
int ndiffs(NODE *s1, NODE *s2)
{
  register int i, diffs=0.0;

  if(s1!=NULL && s2!=NULL)  /* neither == NULL */
  {
    switch(mut_model)
    {
    case infinite:
      i = s1->mutations - s2->mutations;
      return( i>=0 ? i : -i );
    case finite_gamma:
    case finite_flat:
      for(i=0; i<n_sites; i++)
	if(s1->state.seq[i]!=s2->state.seq[i])
	  diffs++;
      return(diffs);
    case stepwise:
      i = s1->state.steps - s2->state.steps;
      return( i>=0 ? i : -i );
    }
  }
  /* if we get this far then at least one is NULL */
  if(s2!=NULL)   /* s1 == NULL; s2 != NULL */
  {
    switch(mut_model)
    {
    case infinite:
      return(s2->mutations);
    case finite_gamma:
    case finite_flat:
      for(i=0; i<n_sites; i++)
	if(s2->state.seq[i] != INIT_SITEVAL)
	  diffs++;
      return(diffs);
    case stepwise:
      return(s2->state.steps);
    }
  }
  if(s1!=NULL)   /* s1 != NULL; s2==NULL */
  {
    switch(mut_model)
    {
    case infinite:
      return(s1->mutations);
    case finite_gamma:
    case finite_flat:
      for(i=0; i<n_sites; i++)
	if(s1->state.seq[i] != INIT_SITEVAL)
	  diffs++;
      return(diffs);
    case stepwise:
      return(s1->state.steps);
    }
  }             
  /* both are NULL */
    return(0);
}

/* count segregating sites */
int getsegregating(void)
{
  int i, segregating=0;
  
  switch(mut_model)
  {
  case stepwise:
    break;
  case infinite:
    segregating = n_sites;
    break;
  case finite_gamma:
  case finite_flat:
    for(i=0; i<n_sites; i++)
      if(hit[i] > 0)
	++segregating;
    break;
  }
  return(segregating);
}

#define KEEPSIZE 100
/* log of factorial */
double lnfact(int n)
{
  static int first_call = 1;
  static double a[KEEPSIZE];
  int j;

  if(first_call)
  {
    first_call = 0;
    for(j=0; j<KEEPSIZE; j++)
      a[j] = 0.0;
  }
#ifndef NDEBUG  
  assert(n >= 0);
#endif

  if(n <= 1)
    return(0.0);
  if (n >= KEEPSIZE)
    return(gammln(n+1.0));
  if(a[n] == 0.0)
    a[n] = gammln(n + 1.0);
  return a[n];
}

/* poisson probability function */
double poisson(int x, double mean)
{
  double lnprob;

#ifndef NDEBUG  
  assert(x >= 0);
  assert(mean >= 0.0);
#endif
  switch(x)   /* for small x, do it the simple way */
  {
  case 0:
    return(exp(-mean));
  case 1:
    return(mean*exp(-mean));
  case 2:
    return(mean*mean*exp(-mean)/2.0);
  case 3:
    return(mean*mean*mean*exp(-mean)/6.0);
  }

  /* for larger x, work with logs to avoid large numbers */
  lnprob = x*log(mean) - mean - lnfact(x);
  return(exp(lnprob));
}

/* binomial coefficient */
double choose(int n, int k)
{
#ifndef NDEBUG  
  assert(n >= 0);
  assert(k >= 0);
  assert(k <= n);
#endif
  return(floor( 0.5 + exp(lnfact(n) - lnfact(k) - lnfact(n-k)) ));
}

/* log gamma, this came from Numerical Recipes */
double   gammln(double xx)
{
  double  x, tmp, ser;
  static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
  -1.231739516, 0.120858003e-2, -0.536382e-5};
  int     j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j = 0; j <= 5; j++)
  {
    x += 1.0;
    ser += cof[j] / x;
  }
  return -tmp + log(2.50662827465 * ser);
}
/* Poisson random deviate, from Numerical Recipes */
double   poidev(double xm)
{
  static double sq, alxm, g, oldm = (-1.0);
  double   em, t, y;

  if (xm < 12.0)
  {
    if (xm != oldm)
    {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1;
    t = 1.0;
    do
    {
      em += 1.0;
      t *= uni();
    } while (t > g);
  } else
  {
    if (xm != oldm)
    {
      oldm = xm;
      sq = sqrt(2.0 * xm);
      alxm = log(xm);
      g = xm * alxm - gammln(xm + 1.0);
    }
    do
    {
      do
      {
	y = tan(Pi * uni());
	em = sq * y + xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
    } while (uni() > t);
  }
  return (em);
}

/****************************************************************
Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small 
a can be without causing trouble.  f(x) doesn't seem to have a series
expansion around x=0.
****************************************************************/
double gamma_dev(double a)
{
  int ia;
  double u, b, p, x, y=0.0, recip_a;

  if(a <= 0)
  {
    fprintf(stderr,"\ngamma_dev: parameter must be positive\n");
    exit(1);
  }

  ia = floor(a);  /* integer part */
  a -= ia;        /* fractional part */
  if(ia > 0)
  {
    y = igamma_dev(ia);  /* gamma deviate w/ integer argument ia */
    if(a==0.0)
      return(y);
  }

  /* get gamma deviate with fractional argument "a" */
  b = (M_E + a)/M_E;
  recip_a = 1.0/a;
  for(;;)
  {
    u = uni();
    p = b*u;
    if(p > 1)
    {
      x = -log( (b-p)/a );
      if( uni() > pow(x, a-1))
	continue;
      break;
    }else
    {
      x = pow(p, recip_a);
      if( uni() > exp(-x))
	continue;
      break;
    }
  }
  return(x+y);
}

/****************************************************************
gamma deviate for integer shape argument.  Code modified from pp
292-293 of Numerical Recipes in C, 2nd edition.
****************************************************************/
double igamma_dev(int ia)
{
  int j;
  double am,e,s,v1,v2,x,y;
  
  if (ia < 1)
  {
    fprintf(stderr,"\nerror: arg of igamma_dev was <1\n");
    exit(1);
  }
  if (ia < 6)
  {
    x=1.0;
    for (j=0; j<ia; j++)
      x *= uni();
    x = -log(x);
  }else
  {
    do
    {
      do
      {
	do
	{                         /* next 4 lines are equivalent */
	  v1=2.0*uni()-1.0;       /* to y = tan(Pi * uni()).     */
	  v2=2.0*uni()-1.0;
	}while (v1*v1+v2*v2 > 1.0);
	y=v2/v1;
	am=ia-1;
	s=sqrt(2.0*am+1.0);
	x=s*y+am;
      }while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    }while (uni() > e);
  }
  return(x);
}

/* is the integer even? */
int is_even(int i)
{
  return( 2*(i/2) == i );
}

/**************************************************************** 
Sum step variable across all nodes in tree.  On return, m.f[i] = sum
of steps and function returns number of terms in sum.
****************************************************************/
int treesumsteps(int *sum, NODE *node)
{
  int n;

  if(node==NULL)
    return(0);

  n = treesumsteps(sum, node->left) + treesumsteps(sum, node->right);

  if(n > 0)       /* We're not at tip. */
    return(n);    /* Only tip nodes count. */

  *sum += node->state.steps;

  return(1);
}

/**************************************************************** 
Calculate moments about the value "origin" by traversing the tree.
On return, m.f[i] = sum of (steps-origin)^i and function returns number
of terms in sum.  Moments are calculated by:

  DblVec *m;
  double origin;

  m = newdblvec(1,4);
  for(i=m->lo; i<=m->hi; i++)  
    m->f[i] = 0.0;
  origin = 0.0       
  n = treemoments(m, origin, tree);
  for(i=m->lo; i<=m->hi; i++)
    m->f[i] /= n;

****************************************************************/
int treemoments(DblVec *m, double origin, NODE *node)
{
  int i, lo=1, n=0;
  double v=1.0, diff;

  if(node==NULL)
    return(0);

  if(node->left != NULL)
    n += treemoments(m, origin, node->left);

  if(node->right != NULL)
    n += treemoments(m, origin, node->right);

  if(n > 0)       /* We're not at tip. */
    return(n);    /* Only tip nodes count. */

  if(m->lo > 1)
    lo = m->lo;
#ifndef NDEBUG
  assert(lo <= m->hi);
#endif
  diff = node->state.steps - origin;

  for(i=1; i < lo; i++)          /* set v = (steps^(lo-1) */
    v *= diff;

  for(i=lo; i <= m->hi; i++)
  {
    /******************************************
    on successive iterations of loop, v equals
    diff, diff^(lo), diff^(lo+1), and so on
    *******************************************/
    v *= diff;
    m->f[i] += (double) v;  /* sum of steps^(i+1) */
  }
  return(1);
}

/****************************************************************
Convert moments (about origin) of the distribution across chromosomes
of repeat counts into central moments of the distribution of 
pairwise differences in repeat counts.  On entry, G->f[i] is the
expectation of steps^i.  On return, 

  G->f[1] = 0
  G->f[2] = E[ (x-y)^2 ]
  G->f[3] = 0
  G->f[4] = E[ (x-y)^4 ]
  G->f[5] = 0
  G->f[6] = E[ (x-y)^6 ]

where x and y are the repeat counts of random individuals in the
populations.
****************************************************************/
void pairwise_moments(DblVec *G)
{
#ifndef NDEBUG  
  assert(G->lo <= 1);
  assert(G->hi >= 2);
#endif

  if(G->hi >= 6)  /* 6th pairwise moment */
  {
    G->f[6] = 2.0 * G->f[6] - 12.0 * G->f[1] * G->f[5]  
              + 30.0 * G->f[2] * G->f[4] 
            - 20.0 * (G->f[3]) * (G->f[3]) ;
  }

  if(G->hi >= 4)  /* 4th pairwise moment: */
  {
    G->f[4] = 6.0 * G->f[2] * G->f[2] 
      - 8.0 * G->f[1] * G->f[3] + 2.0 * G->f[4];
  }

  /* 2 * 2nd moment: */
  G->f[2] = 2.0 * G->f[2] - 2.0 * G->f[1] * G->f[1];

  /* odd moments */
  G->f[5] = G->f[3] = G->f[1] = 0.0;

}

/**************************************************************** 
On input: G->f[2] and G->f[4] contain the 2nd and 4th moments of
the distribution of (x-y) the pairwise differences in step counts.
On return, F[0] and F[1] contain method-of-moments estimates 
of theta0 and tau, assuming the one-step mutational model.
****************************************************************/
void stepwise_mom(double *theta0, double *tau, DblVec *G)
{
#ifndef NDEBUG  
  assert(G->lo <=2);
  assert(G->hi >=4);
#endif

  /* theta0 = sqrt( (G4 - G2)/3 - G2^2 ) */
  *theta0 = ((G->f[4]) - (G->f[2]))/3.0 - (G->f[2]) * (G->f[2]) ;
  if( *theta0 < 0.0 )
    *theta0 = 0.0;           /* avoid imaginary result */
  else
    *theta0 = sqrt( *theta0 );

  /* tau = theta0 - F->f[0] */
  *tau = G->f[2] - *theta0;
}
     
/* convert mutations into steps */
int getsteps(int mutations)
{
  int steps;
  double u, sum, p;

  if(mutations == 0)    /* shortcut */
    return(0);

  u = uni();  /* uniform r.v. */
  sum= 0.0;
  steps = -(mutations+1);
  do{
    sum += p = prob_step(++steps, mutations);
  }while( u > sum && steps < mutations);
  return(steps);
}


/****************************************************************
Return the probability of observing a step value of s given that
x mutations have occurred and each mutational change is either
+1 or -1, the two possibilities have equal probability 1/2.

Under this symetrical 1-step model, the change after one mutation
is 2*b - 1, where b is a Bernoulli random variable.  The change
after x mutations is

  s(x) = 2*B(x) - x

where B(x) is binomial with parameters x and 1/2.  Thus, the
probability that s(x) equals, say, s is equal to the probability that
B(x) equals (x + s)/2.

Note that s(x) + x = 2*B and must therefore be even.  The probability
of s(x) is 0 if this sum is odd.  Also note that 
     
     prob_step(s, x) = prob_step(-s, x)

****************************************************************/
double prob_step(int s, int x)
{
  static int length = 1;
  static double halfpwr[KEEPSIZE] = {1.0};  /* halfpwr[x] = 1/2^x */
  double p2;

#ifndef NDEBUG  
  assert(x >= 0);
#endif
  /* return Pr[s|x] = 0 unless s+x is even */
  if(!is_even(s+x))
    return(0.0);

  if(s < -x || s > x)     /* out of range */
    return(0.0);

  /* store 1/2^x values to minimize calls to pow function */
  while(length <= x && x < KEEPSIZE)
  {
    halfpwr[length] = 0.5 * halfpwr[length-1];
    length += 1;
  }
  if(x < length)
    p2 = halfpwr[x];
  else
    p2 = pow(0.5, (double) x);

  return(choose(x, (x+s)/2) * p2 );
}

/* describe the model of mutation */
void prstepwise(FILE *fp)
{
  fprintf(fp,"\n#Stepwise mutation model with mutational distribution:");
  fprintf(fp,"\n Pr[-1] = Pr[+1] = 1/2");
}
