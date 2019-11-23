/****************************************************************
mmest: Estimate parameters using mismatch distribution
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memstack.h"
#include "mismatch.h"
#include "estimate.h"
#include "bye.h"
#include "header.h"
#define PROGRAM "MMEST"  

/******* prototypes *******/
void usage(void);

/* Main routine */
int    main(int argc, char **argv)
{
  int msize = 0; /* size of mismatch distribution */
  int     i;
  double  eps, epsp1, unity;
  SIMULATION observation, *obs = &observation;
  FILE   *fp, *ifp=stdin;    /* default input is stdin */
  unsigned eflags;            /* flags to indicate what will be estimated */
  unsigned flags;

  /* initialize eflags */
  eflags = THETA0|THETA1|TAU|MSE|ROUGHNESS|SEG ;

  obs->h = NULL;
  obs->smpsiz = 0;
  fp = stdout;

  /* calculate machine epsilon */
  eps = unity = 1.0;
  do{
    eps *= 0.5;
    epsp1 = eps + unity;
  }while(epsp1 > unity);

/****** Command line arguments *********/
  for(i=1; i<argc; i++)
    if(argv[i][0] == '-') switch(argv[i][1])  /* flag args */
    {
    case 'l':
      msize = atoi(argv[i]+2);     /* set length of mismatch distribution */
      break;
    default:
      usage();
    }else         /* non-flag args are input file names */
    {
      ifp = mustopen(argv[i], "r");  /* open input file */
      obs->fname = argv[i];
    }
  if(ifp == stdin)
    obs->fname = "stdin";
  
/****** Print header ***********/
  header(PROGRAM, "(Estimation from Mismatch Distribution)", fp);
  echo_cmdline(argc, argv);

/***** analyze data **********/
  obs->h = NULL;
  obs->e.seg = -1;
  obs->nsites = -1;
  if((msize=input(obs, NULL, NULL, NULL, msize, &flags, ifp, fp)) == 0)
    error("No data");
  if(obs->h == NULL)
    error("Input data must provide histogram");
  if((eflags&SEG) && obs->e.seg == -1)
    error("Input data must provide number of segregating sites");
  sputval("InputFile", 1, &(obs->fname), fp, NO_COMMENT);
  iputval("NSequences", 1, &(obs->smpsiz), fp, 0);
  if(obs->nsites > 0)
    iputval("NSites", 1, &(obs->nsites), fp, NO_COMMENT);
  iputval("Mismatch", msize, obs->h, fp, NO_COMMENT);
  getcumulants(obs->h, msize, obs->e.cumulant); /* calculate cumulants */
  fprintf(fp,"\n\n%c%11s %1s %13s %13s %13s", START_COMMENT, "", "",
	  "mean", "var", "E[(x-mean)^3]");
  fprintf(fp,"\n%12s %1s %13g %13g %13g ;",
      "Cumulants", "=", obs->e.cumulant[0], obs->e.cumulant[1],
	  obs->e.cumulant[2]);

  /****************************************************************
  When 4th arg=NULL, estimate will calculate MSE using estimated
  values of theta0, theta1, and tau.
  ****************************************************************/
  estimate(obs, msize, NULL, obs->h, eflags); 

  /** Print estimates  **/
  print_labels(fp, eflags, !PREPEND_COMMENT);
  print_estimates(obs, fp, eflags, !PREPEND_COMMENT);

  /** Set up ranges **/
  fprintf(fp,
       "\n\n%c mmci ranges are specified by: start-val, end-val, increment",
	  START_COMMENT);
  fprintf(fp,"\n%c Edit the following to control mmci:",
	  START_COMMENT);
  fprintf(fp,"\nRangeLog10Theta0 = %d %d %3.2f ;", 0, 3, 0.5);
  fprintf(fp,"\nRangeGrowth = %d %d %d ;", 0, 3, 1);
  fprintf(fp,"\nRangeTau = %d %d %d ;", 2, 12, 2);

  putc('\n', fp);
  if(fp != stdout)
    fclose(fp);
  return 0;
}

void usage(void)
{
  char buff[30];

  fflush(stdout);
  fprintf(stderr,
	  "\nusage:  mmest [options] inputfile\n where options may include:");
  sprintf(buff,"as in input");
  option("-l<x>", "Set length of mismatch distribution to <x>.", buff);
  putc('\n', stderr);
  exit(1);
}

