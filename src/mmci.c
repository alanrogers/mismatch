/****************************************************************
mmci: Confidence region for mismatch distributions
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "memstack.h"
#include "mismatch.h"
#include "bye.h"
#include "unirand.h"
#include "alloc2d.h"
#include "header.h"
#include "eprintf.h"
/* #include "maxlike.h" */

/******* prototypes *******/
void usage(void);
void set_tables(double *from, double *to, double *by, POPHIST *ph, int msize);
void dummy(void);

/** externals **/
int         iterations=1000;    /* number of iterations to perform */
int         init_n_sites = -1;  /* # of sites in finite sites models */
double      shape = 1.0;        /* shape of gamma distribution */
int         verbose = 1;
int         reset_mutation=0;   /* reset gamma mutation rates each time?*/
extern int  total_iterations;   /* total number of iterations in all */
#if 0
double obs_mae;
#endif

/** defined in iscoales.c **/
extern int *mismatch, ***intermatch;
extern int npops, nwi;              

void usage(void)
{
  char buff[30];
  
  fprintf(stderr,
	  "\nusage:  mmci [options] inputfile\n where options may include:");
  sprintf(buff,"%f", shape);
  option("-g<x>", "Set gamma shape parameter to <x>", buff);
  sprintf(buff,"%d",iterations);
  option("-i<x>", "Set iterations to <x>.", buff);
  sprintf(buff,"as in input");
  option("-l<x>", "Set length of mismatch distribution to <x>.", buff);
  option("-mi", "mutation model = infinite sites", NULL);
  option("-mf", "mutation model = finite sites w/ equal rates", NULL);
  option("-mg", "mutation model = finite sites w/ gamma rates", NULL);
  option("-ms", "mutation model = stepwise", NULL);
  sprintf(buff,"%s", YES(verbose));
  option("-r", "Reset gamma-model mutation rates each time?",
	 YES(reset_mutation));
  option("-v", "Toggle verbose mode", buff);
  putc('\n', stderr);
  exit(1);
}

void    main(int argc, char **argv)
{
  int   i, ival, ntot, count, thisfig, rtnval;
  int   msize = 0; /* size of mismatch distribution */
  int   test;   /* count the tests performed */
  double  eps, epsp1, unity, log10theta0, pval;
  SIMULATION simulation, *sim = &simulation;
  SIMULATION observation, *obs = &observation;
  POPHIST *history=NULL;
  FILE   *fp, *ifp=stdin, *hfp;
  char *hfile = "pophist.ini";
  char *progname = "main";
  double **x=NULL;                /* hold simulated estimates */
  MUTATION_MODEL mut_model = infinite;
  unsigned eflags;    /* which values were present in observed data? */
  unsigned sflags;    /* which will be estimated in simulated data? */
  int       dim_s;    /* number of parameters to be estimated */
  /* for confidence interval */
#if 0
  double from[3], to[3], by[3];   /* indices:0=theta0,1=theta1,2=tau */
#else
  double *from, *to, *by;   /* indices:0=theta0,1=theta1,2=tau */
#endif
  double growth;
  double theta0, theta1, tau;
  double *theory_f;
  ASSIGNMENT *a;

  /* which quantities should be estimated in simulated data sets? */
  sflags = (THETA0 | THETA1 | TAU | MSE);

  assert( (THETA0 & THETA1 & TAU & MSE & MAE & ROUGHNESS & SEG) == 0);

  from = (double *) emalloc(progname, 3 * sizeof(double));
  by = (double *) emalloc(progname, 3 * sizeof(double));
  to = (double *) emalloc(progname, 3 * sizeof(double));


  obs->smpsiz = sim->smpsiz = 0;
  obs->tree =sim->tree = NULL;
  fp = stdout;
  initrand(0);  /* initialize random number generator */
  /* default range of values for log10 theta0 */
  from[0] = -(1.0 + 1.0/3.0);
  to[0]   = 1.0 + 1.0/3.0;
  by[0]   = 1.0/3.0;
  /* default range of values for log10( theta1/theta0 )*/
  from[1] = 0.0;
  to[1]   = 7.0;
  by[1]   = 1.0;
  /* default range of values for tau */
  from[2] = 3.0;
  to[2]   = 8.0;
  by[2]   = 1.0;

  /* calculate machine epsilon */
  eps = unity = 1.0;
  do{
    eps *= 0.5;
    epsp1 = eps + unity;
  }while(epsp1 > unity);

/****** Print header ***********/
  header("MMCI", "(Confidence Region for Mismatch Distributions)", fp);

/****** Command line arguments *********/
  for(i=1; i<argc; i++)
  {
    if(argv[i][0] == '-') switch(argv[i][1])  /* flag args */
    {
    case 'g':
      shape = atof(argv[i]+2);
      break;
    case 'i':
      iterations = atoi(argv[i]+2);
      break;
    case 'l':
      msize = atoi(argv[i]+2);     /* set length of mismatch distribution */
      break;
    case 'm':
      switch(argv[i][2])
      {
      case 'i':
	mut_model = infinite;
	break;
      case 'f':
	mut_model = finite_flat;
	break;
      case 'g':
	mut_model = finite_gamma;
	break;
      case 's':
	mut_model = stepwise;
	break;
      default:
	usage();
      }
      break;
    case 'r':
      reset_mutation = TOGGLE(reset_mutation);
      break;
    case 'v':
      verbose = TOGGLE(verbose);
      break;
    default:
      fprintf(stderr,"\nIllegal command line argument: %s", argv[i]);
      usage();
    }else
    {
      if(ifp!=stdin)
      {
	fprintf(stderr,"\nOnly 1 input file is allowed.");
	usage();
      }
      ifp = mustopen(argv[i], "r");  /* open input file */
      obs->fname = argv[i];
    }
  }
  if(ifp == stdin)
    obs->fname = "stdin";
  
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
  obs->smpsiz = 0;
  obs->h = NULL;

  count = 0;
  eflags = 0;
  while((a=getassignment(ifp)) != NULL)
  {
    count++;
    switch(a->lhs)
    {
    case InputFile:
      obs->fname = dupstring(a->val[0]);
      break;
    case Sampsize:
      obs->smpsiz = atoi(a->val[0]);
      break;
    case NSites:
      obs->nsites = init_n_sites = atoi(a->val[0]);
      break;
    case Histogram:
      if(msize <= 0)   /* histogram size set by input data */
	msize = a->length;
      /* otherwise, histogram size is set by parameter msize */
      obs->h = (int *) emalloc(progname, msize*sizeof(int));
      obs->h[msize-1] = 0;
      fprintf(fp,
      "\n%c Using histogram of length %d. Input histogram had length %d.",
	      START_COMMENT, msize, a->length);
      for(i=0; i<a->length; i++)
      {
	ival = atoi(a->val[i]);  /* round to nearest int */
	/* ensure that input value was an int */
	if(fabs(((double) ival) - atof(a->val[i])) > 0.000001)
	{
	  fflush(stdout);
	  fprintf(stderr,"\nError in input: Non-integer value in histogram");
	  fprintf(stderr,"\n  Illegal value: h[%d] = %s\n",
		  i, a->val[i]);
	  exit(1);
	}
	if(i >= msize)           /* last entry of h is sum of large values */
	  obs->h[msize-1] += ival;
	else
	  obs->h[i] = ival;
      }
      for(i=a->length; i<msize; i++)  /* zero unassigned portion of h */
	obs->h[i] = 0;
      break;
    case Cumulants:
      break;
    case Labels:
      for(i=0; i < a->length; i++)
      {
	if(!strcmp(lowercase(a->val[i]), "theta0")){
	  eflags |= THETA0;
	}
	if(!strcmp(lowercase(a->val[i]), "theta1")){
	  eflags |= THETA1;
	}
	if(!strcmp(lowercase(a->val[i]), "tau")){
	  eflags |= TAU;
	}
	if(!strcmp(lowercase(a->val[i]), "mse")){
	  eflags |= MSE;
	}
	if(!strcmp(lowercase(a->val[i]), "mae")){
	  eflags |= MAE;
	}
	if(!strcmp(lowercase(a->val[i]), "roughness")){
	  eflags |= ROUGHNESS;
	}
	if(!strcmp(lowercase(a->val[i]), "seg")){
	  eflags |= SEG;
	}
      }
      break;
    case Estimates:
      if(eflags==0)
	error("Input file must specify labels before estimates");
      i=0;
      if(eflags & THETA0)
	obs->e.theta0 = atof(a->val[i++]);
      if(eflags & THETA1)
	obs->e.theta1 = atof(a->val[i++]);
      if(eflags & TAU)
	obs->e.tau = atof(a->val[i++]);
      if(eflags & MSE)
	obs->e.mse = atof(a->val[i++]);
      if(eflags & MAE)
	obs->e.mae = atof(a->val[i++]);
      if(eflags & ROUGHNESS)
	obs->e.roughness = atof(a->val[i++]);
      if(eflags & SEG)
	obs->e.seg = atoi(a->val[i++]);
      assert(i == a->length);
      break;
    case RangeLog10Theta0:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[0] = atof(a->val[0]);
      to[0]   = atof(a->val[1]);
      by[0]   = atof(a->val[2]);
      break;
    case RangeGrowth:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[1] = atof(a->val[0]);
      to[1]   = atof(a->val[1]);
      by[1]   = atof(a->val[2]);
      break;
    case RangeTau:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[2] = atof(a->val[0]);
      to[2]   = atof(a->val[1]);
      by[2]   = atof(a->val[2]);
      break;
    case Eof:
      count--;
      break;
    case Badval:
      count--;
      error("input(): Bad value in switch");
    default:
      count--;
      break;
    }
    freeassignment(a);
  }

  if(count == 0)
    error("No assignments in input file.");
  if(obs->smpsiz == 0)
    error("Input data must specify sampsize");
  if(mut_model != infinite && init_n_sites== -1)
    error(
       "Input data must specify NSites unless infinite sites model is used");
  if(obs->h==NULL)
    error("Input data must specify histogram");
  if(eflags == 0)
    error("Input data must give parameter estimates");

  obs->nsubs = history->K;  /* number of subdivisions */
  obs->subsiz = (int *) emalloc(progname, obs->nsubs * sizeof(int));
  /* make sample size evenly divisible by obs->nsubs */
  obs->smpsiz = (obs->smpsiz / obs->nsubs) * obs->nsubs ;
  for(i = 0; i< obs->nsubs; i++)
    obs->subsiz[i] = obs->smpsiz / obs->nsubs; /* sizes of subdivisions */
  sputval("InputFile", 1, &(obs->fname), fp, PREPEND_COMMENT);
  iputval("Iterations", 1, &iterations, fp, PREPEND_COMMENT);
  fprintf(fp,"\n%c total sampsize=%d subdivisions=%d",
	  START_COMMENT, obs->smpsiz, obs->nsubs);
  fprintf(stdout,"\n%c subdivision sizes:", START_COMMENT);
  for(i=0; i < obs->nsubs; i++)
    printf(" %d", obs->subsiz[i]);
  iputval("Sampsize", 1, &(obs->smpsiz), fp, PREPEND_COMMENT);
  fprintf(fp,"\n%cInput data included:", PREPEND_COMMENT);

  /* dim_s counts the number of parameters to be estimated */
  dim_s = countbits(sflags);
  fprintf(fp,"\n%csflags has %d bits turned on", START_COMMENT, dim_s);

  set_tables(from, to, by, history, msize);

  switch(mut_model)
  {
  case finite_gamma:
    fprintf(fp,
	    "\n%c Mutation model = finite-gamma. sites=%d shape param=%f",
	    START_COMMENT, init_n_sites, shape);
    fprintf(fp,"\n%c Mutation rates will be set %s.",
	    START_COMMENT,
	    (reset_mutation ? "independently for each tree" :
	                      "once at the top"));
    break;
  case finite_flat:
    fprintf(fp,"\n%c Mutation model = finite-flat. sites=%d",
	    START_COMMENT, init_n_sites);
    break;
  case infinite:
    fprintf(fp,"\n%c Mutation model = infinite sites", START_COMMENT);
    break;
  default:
    fprintf(stderr,"\nILLEGAL MUTATION MODEL. (Value = %d)", (int) mut_model);
    exit(1);
  }
  init_mutation(mut_model, init_n_sites, shape, reset_mutation);
  fflush(fp);

  /**** allocate *******/
  x = (double **) alloc2d(dim_s, iterations, sizeof(double));
  if(x==NULL)
    error("in main(): Can't allocate x");
  sim->h = (int *) emalloc(progname, msize * sizeof(int));
  theory_f = (double *) emalloc(progname, msize * sizeof(double));


  /* initialize sim */
  sim->tree = NULL;
  sim->smpsiz = obs->smpsiz;
  sim->nsubs = obs->nsubs;
  sim->subsiz = (int *) emalloc(progname, sim->nsubs * sizeof(int));
  for(i = 0; i< obs->nsubs; i++)
    sim->subsiz[i] = obs->subsiz[i];
# if 0
  /* is one of the following needed? */
  set_parameters(&(sim->p), &(obs->p)); /* set sim->p = obs->p */
  set_estimates(&(sim->e), &(obs->e));  /* set sim->e = obs->e */
#endif
      
  /**** set memory stack ****/
  i = (2*sim->smpsiz -1)*sizeof(NODE);
  if(mut_model != infinite)
    i += (2*sim->smpsiz -1)*init_n_sites*sizeof(STATE);
  i += 4000;
  if( setmemstack( i )==NULL )
    error("Can't set memstack");

  /** Print observations to be used in test  **/
  fprintf(fp,"\n%cInput data values:", START_COMMENT);
  print_labels(fp, eflags, PREPEND_COMMENT);
  print_estimates(obs, fp, eflags, PREPEND_COMMENT);
  fprintf(fp,"\n%cTests use following observed values:", START_COMMENT);
  print_labels(fp, sflags, PREPEND_COMMENT);
  print_estimates(obs, fp, sflags, PREPEND_COMMENT);

  /* make sure that eflags contains everything in sflags */
  if( ((eflags ^ sflags) & sflags) > 0)
    error("main: Some bit is turned on in sflags but not in eflags");


  /* print ranges in format required by ci2ptx */
  fprintf(fp,"\nRangeLog10Theta0 = %f %f %f ;",
          from[0], to[0], by[0]);
  fprintf(fp,"\nRangeGrowth = %f %f %f ;",
          from[1], to[1], by[1]);
  fprintf(fp,"\nRangeTau = %f %f %f ;",
          from[2], to[2], by[2]);

  if(verbose)
  {
    /***** Count total number of iterations to be performed ******/
    total_iterations=0;
    for(growth = from[1]; growth < to[1]+0.5*by[1] ; growth += by[1])
      for(tau = from[2]; tau < to[2]+0.5*by[2]; tau += by[2]) /* vary tau */
      {
	for(log10theta0=from[0]; log10theta0 < to[0]+0.5*by[0];
	    log10theta0+= by[0])
	  ++total_iterations;
	/* multiple values of tau unnecessary if growth = 0 */ 
	if(growth == 0.0)
	  break;
      }
    total_iterations *= iterations;
    fprintf(stderr,"\nDoing %d iterations in all\n", total_iterations);
  }


  /***** Confidence interval *********/
  bold_comment("Confidence Region", fp, START_COMMENT);

  thisfig = ntot = 0;
  test = 0;
  for(growth = from[1]; growth < to[1]+0.5*by[1] ; growth += by[1]) 
  {                                                   /* vary growth */
    fprintf(fp,"\n\n%c Begin figure %d", START_COMMENT, ++thisfig);
    fprintf(fp,"\n%10s = %13.5f ;", "Growth", growth);
      fprintf(fp,"\n%c %9s   %13s %13s %13s", START_COMMENT, "", "tau",
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

	/* f_hist gets the theoretical mismatch distribution.
         * fixsum adjusts the final entry so that theory_f sums
         * to unity.
         */
	(void) fixsum(f_hist(history, theory_f, msize), msize);
	if(sflags &MSE) {
	  /*
	   * Here is where obs->e.mse gets clobbered.
	   */
	  fprintf(fp,"\n%c before assignment:obs->e.mse = %g",
		  START_COMMENT, obs->e.mse);
	  obs->e.mse = get_mse(theory_f, obs->h, msize);
	  fprintf(fp,"\n%c after assignment:obs->e.mse = %g",
		  START_COMMENT, obs->e.mse);
	}
	if(sflags & MAE)
	  obs->e.mae = get_mae(theory_f, obs->h, msize);

	fprintf(fp,"\n%c before multisim:sim->e.mse = %g",
	  START_COMMENT, sim->e.mse);
	ntot = multisim(sim, history, iterations, x, theory_f, msize,
			mut_model, sflags);
	fprintf(fp,"\n%c after multisim:sim->e.mse = %g",
	  START_COMMENT, sim->e.mse);

	/*** test of hypothesis ***/
	rtnval = (double) testH0(obs, x, iterations, eps, fp, msize, dim_s,
				 eflags, sflags);
	if(rtnval >= 0)
	{
	  pval = (double) rtnval / (double) ntot;
	  fprintf(fp,"\n%10s = %13.7f %13.7f %13.7f ;",
		  "Test", tau, log10theta0,
		  pval);  /* = p-value */
	  /*record_test(tau, log10theta0, growth, pval);*/
	}else  /* error */
	  fprintf(fp,"\n%10s = %13.7f %13.7f %13d ;",
		  "Test", tau, log10theta0, rtnval);
	fflush(fp);
      }
      /* multiple values of tau unnecessary if growth = 0 */ 
      if(growth == 0.0)
	break;
    }
  }
  bold_comment("end of confidence region", fp, START_COMMENT);
  /*prmaxlike(fp);*/
  putc('\n', fp);
  if(fp != stdout)
    fclose(fp);
  putc('\n', stderr);
  exit(0);
}

/****************************************************************
This procedure isolates some messy code that figures out all the
values of theta and of tau that will ever be needed and then calls
init_theory to initialize the external tables in identity.c.  On
input, from[0] is the initial value of theta0, from[1] that of theta1,
and from[2] that of tau.  to[] and by[] are defined similarly.  Returns the
total number of parameter value combinations to be simulated.
****************************************************************/
void set_tables(double *from, double *to, double *by, POPHIST *ph, int msize)
{
  double growth, ltheta0, theta0, theta1, tau;
  double *tau_val, *theta_val, *v, x;
  int  i, j, n_tau, n_theta, n;
  int ntheta0, ngrowth, nhist;
  POPHIST *ph2;
  char *progname = "set_tables";

  fflush(stdout);
  /* 1st step: calculate dimensions */
  ntheta0 = floor((to[0] + 0.5*by[0] - from[0])/by[0]) + 1;
  ngrowth = floor( (to[1] + 0.5*by[1] - from[1])/by[1]) + 1;
  for(nhist=0, ph2=ph; ph2 != NULL; ph2 = ph2->next)
    nhist += 1;
  n = ntheta0 + ntheta0*ngrowth + nhist;

  /* 2nd step: allocate array of theta values */
  v = (double *) emalloc(progname, n * sizeof(double));
  
  /* 3rd step: put values in array */
  i = 0;
  for(ltheta0=from[0];
      ltheta0 < to[0] + 0.5*by[0];
      ltheta0 += by[0])
  {
    theta0 = pow(10.0, ltheta0);
    x = theta0/(double)ph->next->K;
    if(i>=n)
    {
      fprintf(stderr,"\nset_tables 1: attempting to access v[%d]", i);
      fprintf(stderr,"\n              max index is %d", n-1);
      fprintf(stderr,"\n  This is a symptom of a bug in gcc.");
      fprintf(stderr,"\n  To fix it, don't use -O2 or -O3 when");
      fprintf(stderr,"\n  compiling\n");
      exit(1);
    }
    v[i++] = x * ph->next->K;
    for(growth = from[1];
	growth < to[1] + 0.5*by[1] ;
	growth += by[1]) 
    {
      theta1 = theta0*pow(10.0, growth);
      x = theta1 / (double) ph->K;
      if(i>=n)
      {
	fprintf(stderr,"\nset_tables 2: attempting to access v[%d]", i);
	fprintf(stderr,"\n              max index is %d", n-1);
	fprintf(stderr,"\n  This is a symptom of a bug in gcc.");
	fprintf(stderr,"\n  To fix it, don't use -O2 or -O3 when");
	fprintf(stderr,"\n  compiling\n");
	exit(1);
      }
      v[i++] = x * ph->K;
    }
  }
  for(ph2 = ph; ph2 != NULL; ph2 = ph2->next)
    v[i++] = ph2->theta * ph2->K;
  if(i!=n)
  {
    fprintf(stderr,"\nset_tables: Put %d values into array v.", i);
    fprintf(stderr,"\n  Should have put %d values", n);
    fprintf(stderr,"\n  This is a symptom of a bug in gcc.");
    fprintf(stderr,
	    "\n  To fix it, don't use -O or -O2 or -O3 when compiling\n");
    exit(1);
  }

  /* 4th step: sort array */
  qsort((void *) v, (size_t) n, (size_t) sizeof(double),
	(int (*)(const void*,const void*)) compar );
  /* 5th step: count unique values, set n_theta */
  n_theta = n;
  for(i=1; i<n; i++)
    if(v[i] == v[i-1])
      --n_theta;
  /* 6th step: allocate theta_val */
  theta_val = (double *) emalloc(progname, n_theta * sizeof(double));
  /* 7th step: put unique values into array */
  if(n_theta > 0)
    theta_val[0] = v[0];
  for(i=j=1; i<n; i++)
  {
#ifndef NDEBUG
    assert(j < n_theta);
#endif
    if(v[i] != v[i-1])
      theta_val[j++] = v[i];
  }
#ifndef NDEBUG
  assert(j == n_theta);
#endif
  /* final step: */
  free(v);

  /* set n_tau and tau_val */
  n = 0;
  for(tau = from[2]; tau < to[2]+0.5*by[2]; tau += by[2]) /* vary tau */
    ++n;
  for(ph2 = ph; ph2 != NULL && ph2->next != NULL; ph2 = ph2->next)
    ++n;
  v = (double *) emalloc(progname, n * sizeof(double));
  i=0;
  for(tau = from[2]; tau < to[2]+0.5*by[2]; tau += by[2]) /* vary tau */
    v[i++] = tau;
  for(ph2 = ph; ph2 != NULL && ph2->next != NULL; ph2 = ph2->next)
    v[i++] = ph2->tau;
  qsort((void *) v, (size_t) n, (size_t) sizeof(double),
	(int (*)(const void*,const void*)) compar);
  n_tau = n;
  for(i=1; i<n; i++)
    if(v[i] == v[i-1])
      --n_tau;
  tau_val = (double *) emalloc(progname, n_tau * sizeof(double));
  /* 7th step: put unique values into array */
  if(n_tau > 0)
    tau_val[0] = v[0];
  for(i=j=1; i<n; i++)
  {
#ifndef NDEBUG
    assert(j < n_tau);
#endif
    if(v[i] != v[i-1])
      tau_val[j++] = v[i];
  }
#ifndef NDEBUG
  assert(j == n_tau);
#endif
  /* final step: */
  free(v);


  fflush(stdout);
  fprintf(stderr,"\n%c tau values:", START_COMMENT);
  for(i=0; i<n_tau; i++)
    fprintf(stderr," %f", tau_val[i]);

  fprintf(stderr,"\n%c theta values:", START_COMMENT);
  for(i=0; i<n_theta; i++)
    fprintf(stderr," %g", theta_val[i]);
  
  /* now initialize the tables in identity.c */
  init_theory(msize, tau_val, n_tau, theta_val, n_theta);

  free(theta_val);
  free(tau_val);
}

