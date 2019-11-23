/****************************************************************
ci2txt: Make a character plot of a confidence region.
$Id: ci2txt.c 604 2003-10-27 05:31:55Z rogers $
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memstack.h"
#include "mismatch.h"
#include "alloc2d.h"
#include "maxlike.h"
#include "bye.h"
#include "header.h"
#define PROGRAM "CI2TXT"  
#define HLBLWIDTH 8
#define HLBLPRECISION 2
#define VLBLWIDTH 2
#define VLBLPRECISION 2

/******* types ************/
struct charmat {
  int rowdim, coldim;
  char **base;
};

/******* prototypes *******/
void usage(char *message);
char charval(double pval);
struct charmat *cm_alloc(unsigned rowdim, unsigned coldim);
void cm_check_rowdim(struct charmat *m, unsigned row, char *msg);
void cm_check_coldim(struct charmat *m, unsigned col, char *msg);
void cm_print(struct charmat *m);
void cm_putchar(struct charmat *m, unsigned row, unsigned col, unsigned c);
void cm_hline(struct charmat *m, unsigned row, unsigned col, unsigned length);
void cm_putmap(struct charmat *m, unsigned toprow, unsigned leftcol,
	       char **submat, unsigned rowdim, unsigned coldim);
void cm_vfprint(struct charmat *m, unsigned row, unsigned col,
		char *fmt, double value);
void cm_vsprint(struct charmat *m, unsigned row, unsigned col,
		char *fmt, char *svalue);
void cm_hsprint(struct charmat *m, unsigned row, unsigned col, char *s);
void cm_erase(struct charmat *m);
/* make a new charmat */
struct charmat *cm_alloc(unsigned rowdim, unsigned coldim)
{
  struct charmat *m;
  
  m = (struct charmat *) malloc( sizeof(struct charmat) );
  if(m == NULL)
    error("newcharmat: malloc");
  m->rowdim = rowdim;
  m->coldim = coldim;
  m->base = (char **) alloc2d(m->rowdim, m->coldim, sizeof(char));
  if(m->base == NULL)
    error("newcharmat: alloc2d");
  return(m);
}

void cm_check_rowdim(struct charmat *m, unsigned row, char *msg)
{
  if(row >= m->rowdim)
  {
    fprintf(stderr,"\nError: %s: exceeded rowdim", msg);
    fprintf(stderr,"\nrow=%d rowdim=%d\n", row, m->rowdim);
    exit(1);
  }
}

void cm_check_coldim(struct charmat *m, unsigned col, char *msg)
{
  if(col >= m->coldim)
  {
    fprintf(stderr,"\nError: %s: exceeded coldim", msg);
    fprintf(stderr, "\ncol=%d coldim=%d\n", col, m->coldim);
    exit(1);
  }
}

/* print plot */
void cm_print(struct charmat *m)
{
  int i, j;

  for(i=0; i<m->rowdim; i++)
  {
    putchar('\n');
    for(j=0; j<m->coldim; j++)
      putchar(m->base[i][j]);
  }
}

/* put a character onto a plot */
void cm_putchar(struct charmat *m, unsigned row, unsigned col, unsigned c)
{
  cm_check_rowdim(m, row, "cm_putchar");
  cm_check_coldim(m, col, "cm_putchar");
  m->base[row][col] = c;
}

/* draw a horizontal line on a character buffer */
void cm_hline(struct charmat *m, unsigned row, unsigned col, unsigned length)
{
  int i;

  cm_check_rowdim(m, row, "cm_hline");
  cm_check_coldim(m, col+length-1, "cm_hline");
  cm_putchar(m, row, col, '+');
  for(i=1; i<length-1; i++)
    cm_putchar(m, row, col+i, '-');
  cm_putchar(m, row, col+length-1, '+');
}

char charval(double pval)
{
  int ip;

  /* re-express as an integer to avoid rounding errors */
  ip = floor(1000.0*pval + 0.5);
  if(ip < 10)
    return('O');  /* reject at 1% level */
  if(ip < 50)
    return('.');  /* w/i 99% confidence region */
  return('*');   /* w/i 95% confidence region */
}

/* copy rectangular character map onto a character matrix */
void cm_putmap(struct charmat *m, unsigned toprow, unsigned leftcol,
	       char **submat, unsigned rowdim, unsigned coldim)
{
  int i, j;

  cm_check_rowdim(m, toprow + rowdim + 1, "cm_putmap");
  cm_check_coldim(m, leftcol + 2*coldim, "cm_putmap");
  cm_hline(m, toprow, leftcol, 2*coldim+1); 
  for(i=0; i<rowdim; i++)
  {
    cm_putchar(m, toprow+i+1, leftcol, '|');
    for(j=0; j<coldim; j++)
      cm_putchar(m, toprow+i+1, leftcol+1+2*j,
		 (unsigned) submat[rowdim-i-1][j]);
    cm_putchar(m, toprow+i+1, leftcol+1+2*coldim-1, '|');
  }
  cm_hline(m, toprow+rowdim+1, leftcol, 2*coldim+1);
}

/* print double value vertically */
void cm_vfprint(struct charmat *m, unsigned row, unsigned col,
		char *fmt, double value)
{
  char svalue[20];
  int i, len;

  cm_check_coldim(m, col, "cm_vfprint");
  sprintf(svalue,fmt, value);
  len = strlen(svalue);
  cm_check_rowdim(m, row + len - 1, "cm_vfprint");
  for(i=0; i<len; i++)
    cm_putchar(m, row+i, col, (unsigned) svalue[i]);
}

/* print string vertically */
void cm_vsprint(struct charmat *m, unsigned row, unsigned col,
		char *fmt, char *svalue)
{
  int i, len;

  len = strlen(svalue);
  cm_check_coldim(m, col, "cm_vsprint");
  cm_check_rowdim(m, row+len-1, "cm_vsprint");

  for(i=0; i<len; i++)
    cm_putchar(m, row+i, col, (unsigned) svalue[i]);
}

/* print a string horizontally */
void cm_hsprint(struct charmat *m, unsigned row, unsigned col, char *s)
{
  int i, len;

  len = strlen(s);
  cm_check_rowdim(m, row, "cm_hsprint");
  cm_check_coldim(m, col+len-1, "cm_hsprint");
  for(i=0; i<strlen(s); i++)
    cm_putchar(m, row, col+i, s[i]);
}

/* erase a character matrix */
void cm_erase(struct charmat *m)
{
  int i, j;

  for(i=0; i<m->rowdim; i++)
    for(j=0; j<m->coldim; j++)
      cm_putchar(m, i, j, ' ');
}

void usage(char *message)
{
  fprintf(stderr,"\nError: %s\n", message);
  fprintf(stderr,"\nUsage: ci2txt inputfile\n");
  exit(1);
}

/* Main routine */
void    main(int argc, char **argv)
{
  int     i, plotwidth, plotheight;
  int     n_theta0, n_growth, n_tau;
  int     igrowth, itheta0, itau;
  int     hlblwidth=HLBLWIDTH, hlblprecision=HLBLPRECISION;
  int     vlblwidth=VLBLWIDTH, vlblprecision=VLBLPRECISION;
  FILE   *ifp=stdin;
  char    hlblfmt[20], vlblfmt[20];
  double  growth_from, growth_to, growth_by, need_growth;
  double  l10theta0_from, l10theta0_to, l10theta0_by, need_l10theta0;
  double  tau_from, tau_to, tau_by, need_tau;
  double  log10theta0, growth, tau, pval;
  char ***map, buff[20];
  struct charmat *m;

  ASSIGNMENT *a;

/****** Print header ***********/
  header(PROGRAM, "(Text display of confidence region)", stdout);

  /** command line arguments **/
  for(i=1; i<argc; i++)
  {
    if(argv[i][0] == '-')
      usage(argv[i]);
    else
    {
      printf("\nInputFile = %s", argv[i]);
      if(ifp!=stdin)
	usage("Too many input files");
      ifp = (FILE *) mustopen(argv[i], "r");  /* open input file */
    }
  }

  sprintf(hlblfmt, "%%%d.%df", hlblwidth, hlblprecision);
  sprintf(vlblfmt, "%%%d.%dg", vlblwidth, vlblprecision);
  /***** read ranges from data **********/
  need_l10theta0 = need_growth = need_tau = 1;
  do{
    a=getassignment(ifp);
    if(a==NULL)
      error("Ran out of data before finding ranges");
    switch(a->lhs)
    {
    case RangeLog10Theta0:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      l10theta0_from = atof(a->val[0]);
      l10theta0_to   = atof(a->val[1]);
      l10theta0_by   = atof(a->val[2]);
      need_l10theta0 = 0;
      break;
    case RangeGrowth:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      growth_from = atof(a->val[0]);
      growth_to   = atof(a->val[1]);
      growth_by   = atof(a->val[2]);
      need_growth = 0;
      break;
    case RangeTau:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      tau_from = atof(a->val[0]);
      tau_to   = atof(a->val[1]);
      tau_by   = atof(a->val[2]);
      need_tau = 0;
      break;
    default:
      error("Illegal input while reading ranges");
    }
    freeassignment(a);
  }while(need_growth || need_l10theta0 || need_tau);

  if(need_l10theta0)
    fprintf(stderr,"\nNeed RangeLog10Theta0");
  if(need_growth)
    fprintf(stderr,"\nNeed RangeGrowth");
  if(need_tau)
    fprintf(stderr,"\nNeed RangeTau");
  if(need_l10theta0 || need_growth || need_tau)
    error("Must specify ranges in input file");

  printf("\nRangeLog10Theta0: %f %f %f ",
	  l10theta0_from, l10theta0_to, l10theta0_by);
  printf("\nRangeGrowth: %f %f %f ",
	  growth_from, growth_to, growth_by);
  printf("\nRangeTau: %f %f %f ",
	  tau_from, tau_to, tau_by);

  /* find dimensions */
  n_theta0 = (l10theta0_to - l10theta0_from)/l10theta0_by + 1.5;
  n_growth = (growth_to - growth_from)/growth_by + 1.5;
  n_tau = (tau_to - tau_from)/tau_by + 1.5;
  plotwidth = 2*n_tau+2+hlblwidth+2;
  plotheight = n_theta0+vlblwidth+5;
  printf("\nplotheight: %d  plotwidth: %d", plotheight, plotwidth);

  /* allocate character map */
  map = (char ***) mustalloc(n_growth * sizeof(char **));

  /* initialize map */
  for(igrowth=0; igrowth<n_growth; igrowth++)
  {
    map[igrowth] = (char **) alloc2d(n_theta0, n_tau+1, sizeof(char));
    if(map[igrowth] == NULL)
      error("alloc2d");
    for(itheta0=0; itheta0<n_theta0; itheta0++)
    {
      for(itau=0; itau<n_tau; itau++)
	map[igrowth][itheta0][itau] = ' ';
      map[igrowth][itheta0][n_tau] = '\0';
    }
  }

  /* allocate plot */
  m = (struct charmat*) cm_alloc(plotheight, plotwidth);

  printf("\nVarying log10theta0 from %g to %g by %g (%d values)",
	  l10theta0_from, l10theta0_to, l10theta0_by, n_theta0);
  printf("\nVarying log10(theta1/theta0) from %g to %g by %g (%d values)",
	  growth_from, growth_to, growth_by, n_growth);
  printf("\nVarying tau from %g to %g by %g (%d values)",
	  tau_from, tau_to, tau_by, n_tau);

  /***** Confidence interval *********/
  printf("\nDefinitions:");
  printf("\n %c : Reject at 1 percent significance level", charval(0.0));
  printf("\n %c : Within 99 percent confidence region", charval(0.04));
  printf("\n %c : Within 95 percent confidence region", charval(0.1));
  while( ( a = getassignment(ifp) ) != NULL)
  {
    if(a==NULL)
      break;
    switch(a->lhs)
    {
    case Growth:  /* new value of growth: start new plot */
      growth = atof(a->val[0]);
      igrowth = floor((growth - growth_from)/growth_by + 0.5);
      break;
    case Test:  /* print a test result */
      if(a->length != 3)
      {
	fflush(stdout);
	fprintf(stderr,
		"\nInput error: Each test should be a vector of 3 vals.");
	fprintf(stderr, "\n Bad line: Test = ");
	if(a->length > 0)
	  fprintf(stderr,"%s", a->val[0]);
	for(i=1; i<a->length; i++)
	  fprintf(stderr," %s", a->val[i]);
	fprintf(stderr," ;\n");
	exit(1);
      }
      tau = atof(a->val[0]);
      if(growth == 0)
	tau = tau_from + 0.5*(tau_to - tau_from);
      itau = floor((tau - tau_from)/tau_by + 0.5);
      log10theta0 = atof(a->val[1]);
      itheta0 = floor((log10theta0 - l10theta0_from)/l10theta0_by + 0.5);
      pval = atof(a->val[2]);
      if(pval < 0.0)
      {
	printf("\n%%     %11s   %8.5g %13.5g    p-val=%g",
		"ERROR",  tau, log10theta0, pval);
	continue;
      }
      map[igrowth][itheta0][itau] = charval(pval);
      break;
    default:
      error("Illegal lhs reading confidence interval");
    }
    freeassignment(a);
  }

  /** print map **/
  for(igrowth=0; igrowth<n_growth; igrowth++)
  {
    cm_erase(m);   /* erase all characters in plot */
    growth = growth_from + igrowth*growth_by;
    printf("\n\nPanel %d: N1/N0 = 10^%g:\n", igrowth, growth);
    cm_putmap(m, 0, hlblwidth+1+2, map[igrowth], n_theta0, n_tau);
    cm_vsprint(m, 0, 0, "%s", "log theta0");
    for(itheta0=0; itheta0<n_theta0; itheta0++)
    {
      sprintf(buff, hlblfmt, l10theta0_from + itheta0*l10theta0_by);
      cm_hsprint(m, n_theta0 - itheta0, 1, buff);
    }
    if(growth != 0.0)
    {
      for(itau=0; itau<n_tau; itau++)
	cm_vfprint(m, n_theta0+2, hlblwidth+4+2*itau, vlblfmt,
		   tau_from + itau*tau_by);

      cm_hsprint(m, n_theta0+3+vlblwidth,
		 hlblwidth+4+2*n_tau-1 - strlen("tau"), "tau");
    }
    else
      cm_hsprint(m, n_theta0+3,
		 hlblwidth+4+2*n_tau-1 - strlen("undefined"),
		 "undefined");
    cm_print(m);
  }
  putchar('\n');
  exit(0);
}

