/****************************************************************
mmcistr: Confidence interval for mismatch distributions from STR data
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "memstack.h"
#include "mismatch.h"
#include "unirand.h"
#include "alloc2d.h"
#include "bye.h"
#include "header.h"
#include "readstr.h"
#include "str_misc.h"
#include "samedist.h"
/* #include "maxlike.h" */
#define MAXGROUPS 10

/******* prototypes *******/
void usage(void);
void unitsum(int *source, double *target, int length);
int from_to_by(double *from, double *to, double *by, char *p1);

/** externals **/
int     simulations=100;     /* number of simulations to perform */
int     randomizations=100;   /* used by samedist */
double    shape = 2.0;        /* shape of gamma distribution */
int     verbose = 1;

/** defined in iscoales.c **/
extern int *mismatch, ***intermatch;  /* set by alloc_match */
extern int npops, nwi;                /* ...... alloc_match */

/* convert an int vector to a double vector */
void unitsum(int *source, double *target, int length)
{
  int i;
  double sum=0.0;

  for(i=0; i<length; i++)
    sum += source[i];

  for(i=0; i<length; i++)
    target[i] = source[i] / sum;
}
/****************************************************************
Get from, to, and by values from a string whose format is
":from:to:by", where from, to, and by are double numbers in format
xxx.yyy.
****************************************************************/
int from_to_by(double *from, double *to, double *by, char *p1)
{
  char *p2, *p3;

  if(*p1 != ':')
    return(-1);
  p1 += 1;
  p2 = strchr(p1, ':');
  if(p2 == NULL)
    return(-1);
  *p2++ = '\0';   /* now p2 points to 2nd string */

  p3 = strchr(p2, ':');
  if(p3 == NULL)
    return(-1);
  *p3++ = '\0';  /* now p3 points to 3rd string */

  *from = atof(p1);
  *to   = atof(p2);
  *by   = atof(p3);
  return(0);
}

void usage(void)
{
  char buff[30];
  
  fprintf(stderr,
  "\nusage:  mmcistr [options] inputfile\n where options may include:");
  sprintf(buff,"%f", shape);
  option("-g<x>", "Set gamma shape parameter to <x>", buff);
  sprintf(buff,"as in input");
  option("-l<x>", "Set length of mismatch distribution to <x>.", buff);
  option("-pw:x:y:z", "vary param w from x to y by z", NULL);
  option("", "Params are 0:theta0, 1:theta1, 2: tau", NULL);
  sprintf(buff, "%d", randomizations);
  option("-r<x>", "Set number of randomizations", buff);
  sprintf(buff,"%d",simulations);
  option("-s<x>", "Set simulations to <x>.", buff);
  option("-v", "Toggle verbose mode", YES(verbose));
  putc('\n', stderr);
  exit(1);
}

/* Main routine */
void    main(int argc, char **argv)
{
  int gid[MAXGROUPS];  /* gid[i] is identifier of i'th group */
  int smp[MAXGROUPS];  /* smp[i] = number of individuals in group i */
  int ngrps, nloci, maxrepeat, sampsize;
  int ***frq;          /* frq[i][j][k] = # w/ count k in locus j, group i */
  int *f, *ivec;
  INTVEC *g;
  int   i, j, k, it, ntot, thisfig, total_simulations;
  int   test;   /* count the tests performed */
  double  log10theta0, pval;
  NODE *tree;
  POPHIST *history=NULL;
  FILE   *fp, *ifp=stdin, *hfp;
  char *hfile = "pophist.ini";
  /* for confidence interval */
  double from[3], to[3], by[3];   /* indices:0=theta0,1=theta1,2=tau */
  double growth;
  double theta0, theta1, tau;
  double  **obsmm;           /* matrix of observed mismatch distributions */
  double  **simmm;           /* matrix of simulated mismatch distributions */
  double mutation;            /* mutation rate */

  fp = stdout;
  initrand(0);  /* initialize random number generator */
  /* default range of values for log10 theta0 */
  from[0] = 0.0;
  to[0]   = 2.0;
  by[0]   = 0.5;
  /* default range of values for log10( theta1/theta0 )*/
  from[1] = 0.0;
  to[1]   = 7.0;
  by[1]   = 1.0;
  /* default range of values for tau */
  from[2] = 0.0;
  to[2]   = 100.0;
  by[2]   = 10.0;

/****** Print header ***********/
  header("mmcistr",
	 "(Confidence Interval for STR Mismatch Distributions)", fp);

/****** Command line arguments *********/
  for(i=1; i<argc; i++)
  {
    if(argv[i][0] == '-') switch(argv[i][1])  /* flag args */
    {
    case 'g':
      shape = atof(argv[i]+2);
      break;
    case 'i':
      simulations = atoi(argv[i]+2);
      break;
    case 'v':
      verbose = !verbose;
      break;
    case 'p':
      switch(argv[i][2])
      {
      case '0':  /* range for theta0 */
	if(from_to_by(from+0, to+0, by+0, argv[i]+3))
	  usage();
	break;
      case '1':  /* range for theta1 */
	if(from_to_by(from+1, to+1, by+1, argv[i]+3))
	  usage();
	break;
      case '2':  /* range for tau */
	if(from_to_by(from+2, to+2, by+2, argv[i]+3))
	  usage();
	break;
      default:
	usage();
      }
      break;
    case 'r':
      randomizations = atof(argv[i]+2);
      break;
    default:
      usage();
    }else
    {
      if(ifp!=stdin)
      {
	fprintf(stderr,"\nOnly 1 input file is allowed.");
	usage();
      }
      ifp = mustopen(argv[i], "r");  /* open input file */
    }
  }
  
  hfp = fopen(hfile, "r");
  if (hfp == NULL)
  {
    fprintf(stderr, "\nWarning: Can't open history file \"%s\".", hfile);
    exit(1);
  }
  fprintf(stderr,"\nReading history from file %s", hfile);
  history = gethistory(hfp);
  /* make sure history list has at least two levels: */
  if(history->next==NULL)
    history->next = newhistory(history->theta, history->mn,
			       history->tau, history->K);
  writehistory(fp, history, START_COMMENT);

/***** read data **********/
  readstr(ifp, gid, smp, MAXGROUPS, &frq, &ngrps, &nloci, &maxrepeat);

  /* get total sample size */
  for(i=sampsize=0; i<ngrps; i++)
    sampsize += smp[i];

  fprintf(fp,
   "\n%c ngroups=%d nloci=%d maxrepeat=%d overall sample size = %d genes",
	  START_COMMENT, ngrps, nloci, maxrepeat, sampsize);
  fprintf(fp, "\n%c Doing %d simulations of each history",
	  START_COMMENT, simulations);
  fflush(fp);
  init_mutation(stepwise, -1, -1.0, 0);

  /* allocate arrays for mismatch distributions */
  obsmm = (double **) alloc2d(nloci, maxrepeat+1, sizeof(double));
  simmm = (double **) alloc2d(simulations, maxrepeat+1, sizeof(double));
  if(obsmm == NULL || simmm == NULL)
    error("not enough memory for obsmm and simmm");

  /* initialize both arrays with zeroes */
  memset( (void *) obsmm[0], 0, sizeof(double)*nloci*(maxrepeat+1) );
  memset( (void *) simmm[0], 0, sizeof(double)*simulations*(maxrepeat+1) );

  /* calculate observed mismatch distributions */
  f = (int *) mustalloc( (maxrepeat+1) * sizeof(int) );
  ivec = (int *) mustalloc((maxrepeat+1) * sizeof(int) );
  for(j=0; j<nloci; j++)
  {
    /* lump group frequency distributions to get f. */
    /* I should add intermatch comparisons here. */
    memset( (void *) f, 0, (maxrepeat+1)*sizeof(int));
    for(i=0; i<ngrps; i++)
      for(k=0; k<=maxrepeat; k++)
	f[k] += frq[i][j][k];

    get_mismatch(ivec, f, f, maxrepeat);
    /* divide vector by its sum ans store in obsmm */
    unitsum(ivec, obsmm[j], maxrepeat+1);
  }
  
  /**** set memory stack ****/
  i = (2*sampsize -1)*sizeof(NODE) + 4000;
  if( setmemstack( i )==NULL )
    error("Can't set memstack");

  /* print ranges in format required by ci2ptx */
  fprintf(fp,"\nRangeLog10Theta0 = %f %f %f ;",
          from[0], to[0], by[0]);
  fprintf(fp,"\nRangeGrowth = %f %f %f ;",
          from[1], to[1], by[1]);
  fprintf(fp,"\nRangeTau = %f %f %f ;",
          from[2], to[2], by[2]);

  /***** Count total number of simulations to be performed ******/
  total_simulations=0;
  for(growth = from[1]; growth < to[1]+0.5*by[1] ; growth += by[1])
    for(tau = from[2]; tau < to[2]+0.5*by[2]; tau += by[2]) /* vary tau */
    {
      for(log10theta0=from[0]; log10theta0 < to[0]+0.5*by[0];
	  log10theta0+= by[0])
	++total_simulations;
      /* multiple values of tau unnecessary if growth = 0 */ 
      if(growth == 0.0)
	break;
    }
  total_simulations *= simulations;
  fprintf(fp,"\n%c Doing %d simulations in all\n",
	  START_COMMENT, total_simulations);

  /***** Confidence interval *********/
  bold_comment("STR Confidence Interval", fp, START_COMMENT);

  thisfig = ntot = 0;
  test = 0;
  /* allocate vector centered about zero of length maxrepeat+1 */
  g = newintvec(-maxrepeat/2, 1+maxrepeat - maxrepeat/2);
  for(growth = from[1]; growth < to[1]+0.5*by[1] ; growth += by[1]) 
  {                                                   /* vary growth */
    fprintf(fp,"\n\n%cBegin figure %d", START_COMMENT, ++thisfig);
    fprintf(fp,"\n%10s = %13.5f ;", "Growth", growth);
      fprintf(fp,"\n%c%9s   %13s %13s %13s", START_COMMENT, "", "tau",
	      "log10[theta0]", "p-val");
    for(tau = from[2]; tau < to[2]+0.5*by[2]; tau += by[2]) /* vary tau */
    {
      for(log10theta0=from[0];
	  log10theta0 < to[0]+0.5*by[0];
	  log10theta0+= by[0])
      {
	theta0 = pow(10.0, log10theta0);
	history->next->theta = theta0 / (double) history->next->K;
        theta1 = theta0 * pow(10.0, growth);
	history->theta = theta1 / (double) history->K;
        history->tau         = tau;

	for(it=0; it<simulations; it++)
	{
	  /**************************************************************** 
          Draw mutation rate from a gamma distribution with mean 1 and
          shape as specified.
	  ****************************************************************/
	  mutation = gamma_dev(shape)/shape;
	  tree = iscoales(ngrps, smp, history, maxrepeat+1, mutation);
	  /* fill g with zeroes */
	  memset( (void *) g->b, 0, sizeof(int)*(maxrepeat+1) );
	  str_getfrq(g, tree);  /* get frequency distribution */
	  get_mismatch(ivec, g->b, g->b, maxrepeat);
	  /* divide vector by its sum ans store in simmm */
	  unitsum(ivec, simmm[it], maxrepeat+1);
	  stackfree();
	}

	pval = samedist(obsmm, nloci, simmm, simulations, maxrepeat+1,
			randomizations, dist);
	fprintf(fp,"\n%10s = %13.7f %13.7f %13.7f ;",
		"Test", tau, log10theta0,
		pval);  /* = p-value */
	fflush(fp);
      }
      /* multiple values of tau unnecessary if growth = 0 */ 
      if(growth == 0.0)
	break;
    }
  }
  bold_comment("end of STR confidence region", fp, START_COMMENT);
  if(fp != stdout)
    fclose(fp);
  exit(0);
}

