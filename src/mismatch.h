#ifndef MISMATCH_H
#define MISMATCH_H
#define DEBUG
#define Max(X,Y) ((X)>(Y) ? (X) : (Y) )
#define Min(X,Y) ((X)<(Y) ? (X) : (Y) )
#define MAXESTIMATES 7
#define BADDBL ((double) -999.9999)
#define MAXIMUM_ESTIMATE 1e7
#define TOGGLE(i) (!(i))
#define YES(i) ((i) ? "YES" : "NO")
#define START_COMMENT '%'
#define PREPEND_COMMENT 1
#define NO_COMMENT 0
/*
 * Definitions for use with eflags.  Do not change these values.
 */
#define THETA0       1
#define THETA1       2
#define TAU          4
#define MSE          8
#define MAE          16
#define ROUGHNESS    32
#define SEG          64
#include "mytypes.h"
#include "iscoales.h"

/* simulation parameters */
struct param {
  double theta0, theta1, tau;
};

/* quantities estimated from data */
struct estimates {
  double theta0;
  double theta1;
  double tau;
  double roughness;
  double mse;
  double mae;
  int seg;
  double cumulant[3];
};

/**
   This structure is inconsistent with the one of the same name in trimseq.h
**/
typedef struct {
  char *fname;   /* name of data file */
  NODE *tree;
  int   nsubs;   /* subdivisions in ultimate epoch */
  int  *subsiz;  /* size of i'th subdivision in ultimate epoch */
  int   smpsiz;  /* sum of subsiz[i] */
  int   nsites;  /* # of sites in data */
  int  *h;       /* mismatch distribution of length msize */
  struct param p;     /* simulation parameters */
  struct estimates e;   /* estimates */
} SIMULATION;


/** defined in lu.c  */
int     lu_factor(double **a, int n, int *rowpvt, int *colpvt);
int     col_lu_solve(double **a, int n, double *b, int *pvt, double *wa);

/** defined in func.c **/
double *ffun(double *x, double *f);
double **Jfun(double *x, double **J);
void datadef(double *xin);
double *initval(double *p, int *h, int max, double *mom);
double scalarfun(double *p);
double *ffun2(double *p, double *f);

/** defined in getdata.c **/
int    input(SIMULATION *s, double *from, double *to, double *by, int msize,
	     unsigned *eflags, FILE *fp, FILE *ofp);
void iputval(char *name, int len, int *vec, FILE *fp, int comment);
void fputval(char *name, int len, double *vec, FILE *fp, int comment);
void sputval(char *name, int len, char **vec, FILE *fp, int comment);
char *getstring(char *str, int maxlen, FILE *fp);
int instring(int c);
ASSIGNMENT *getassignment(FILE *fp);
int getdouble(double *x, FILE *fp);
void freeassignment(ASSIGNMENT *a);

/** defined in misc.c **/
void  echo_cmdline(int argc, char **argv);
char *dupstring(char *s1);
double quantile(double p, double *vec, int size);
int compar(double *x, double *y);
void prmat(double **m, int rows, int cols, FILE *fp);
void option(char *flag, char *msg, char *def);
void complain(char *s, char *e);
void getcumulants(int *h, int msize, double *m);
int veclength(int *h, int dim);
void pictex_rvec(FILE *fp, int len, double *v, char *lbl);
void bold_comment(char *msg, FILE *fp, int comment_char);
void print_estimates(SIMULATION *obs, FILE *fp, unsigned eflags,
		     int docomment);
void print_labels(FILE *fp, unsigned eflags, int docomment);
unsigned countbits(unsigned u);
char *lowercase(char *s);

/** defined in simulate.c **/
int simulate(SIMULATION *s, POPHIST *history, double *theory_f, int msize,
	     MUTATION_MODEL mut_model, unsigned eflags);
int multisim(SIMULATION *sim, POPHIST *history, int iterations, double **x,
	     double *theory_f, int msize, MUTATION_MODEL mut_model,
	     unsigned sflags);
int copy_estimates(double **x, int col, struct estimates *e, unsigned flags);

/** defined in testH0.c **/
int testH0(SIMULATION *obs, double **x, int iterations, double eps,
	   FILE *fp, int msize, int dim_s, unsigned eflags, unsigned sflags);
void pseudomoments(double **x, int dim, int iterations, 
		   double *med, double **c);
double pseudovar(double *v, int len);
void moments(double **x, int dim, int iterations, 
		   double *mean, double **c);

/** defined in lblaxis.h **/
void lblaxis(double *min, double *max, int *ntic, char **ticmark,
	     double *ticvalue, int logscale);
void pictex_ticks(FILE *fp, int nticks, char **tickmark, double *tickvalue);

#endif /* MISMATCH_H */
