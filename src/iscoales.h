#ifndef ISCOALES_H
#define ISCOALES_H
#define Pi 3.14159265358979323846
#define MAX 1000
#define SQR(x)  ((x)*(x))

typedef enum { infinite, finite_gamma, finite_flat, stepwise } MUTATION_MODEL;

typedef struct list
{
  int     n;			/* # of descendants */
  int    *ndx;			/* array of indices of descendants */
}       LIST;
typedef char SITE;

/* state variable's interpretation depends on mutational model */
/* under infinite sites model, state variable = # of mutations */
typedef union {
  SITE *seq;       /* sequence sites : finite sites models */
  int  steps;      /* # of steps: stepwise model */
} STATE;

typedef struct old_node
{
  float   branch;
  int     pop;                  /* which pop is this node currently in?*/
  int     mutations;
  struct list *d;		/* descendants of this node */
  struct node *left, *right;
} ONODE;

typedef struct node
{
  double   branch;
  int     pop;                  /* which pop is this node currently in?*/
  int     pop0;                 /* which pop is did node start in?*/
  STATE   state;                /* state of chromosome */
  int    mutations;
  /*****************************************************************/
  struct list *d;		/* descendants of this node */
  struct node *left, *right, *ancestor;
} NODE;

typedef struct pophist
{
  double mn;              /* m*n, migrants per generation */
  double theta;           /* 2*N*u where N is subpop size */
  double tau;             /* 2*u*time in generations */
                         /* for infinity, set tau=-1 */
  int K;                 /* # of subpopulations */
  struct pophist *next;  /* next pophist parameters */
} POPHIST ;

/*** prototypes defined in iscoales.c ******/
void  **stackalloc2d(unsigned dim1, unsigned dim2, unsigned size);
NODE   *newnode(NODE * n1, NODE * n2);
void   check_s(char *msg, int *s, int K, NODE **node, int S, int SS);
double   getr(double S, double SS, double mn, double theta);
NODE   *iscoales(int init_nsubs, int *init_subsize, POPHIST *history,
		 int init_msize, double mut_rate);
POPHIST *gethistory(FILE *fp);
POPHIST *newhistory(double theta, double mn, double tau, int K);
void freehistory(POPHIST *ph);
POPHIST *duphist(POPHIST *new, POPHIST *old);
void writehistory(FILE *fp, POPHIST *history, int comment);
int icompar(int *x, int *y);
int get_collapse_vector(int K1, int K2, int *v);
void prtree(NODE * node, int indent);
double treedepth(NODE *node);
double treelength(NODE *node);
int treemutations(NODE *node);
int getmatch(NODE *tree);
void visit(NODE *tree, int total, int origin);
void godown(NODE *tree, int sum, int origin);
void clear_externals(void);
double init_mutation(MUTATION_MODEL mut_model, int init_n_sites,
		   double init_shape, int init_reset);
void set_gamma_rates(void);
void clear_mutation(void);
void mutate(NODE *node, STATE *inherited, double mut_rate);
int getsite(void);
SITE *copyseq(SITE *s1);
int ndiffs(NODE *s1, NODE *s2);
int getsegregating(void);
void    crosstab(NODE * node);
LIST   *newlist(int n, int *ndx);
void prlist(LIST *l);
LIST   *mergelist(LIST * l1, LIST * l2);
void  **stackalloclt(int dim, int size);
double lnfact(int n);
double poisson(int x, double mean);
double choose(int n, int k);
double   gammln(double xx);
double   poidev(double xm);
double   gamma_dev(double a);
double   igamma_dev(int ia);
int    is_even(int i);
int    treesumsteps(int *sum, NODE *node);
int    treemoments(DblVec *m, double origin, NODE *node);
void   pairwise_moments(DblVec *G);
void   stepwise_mom(double *theta0, double *tau, DblVec *G);
int    getsteps(int mutations);
double   prob_step(int s, int x);
void   prstepwise(FILE *fp);
#endif /* ISCOALES_H */
