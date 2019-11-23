/****************************************************************
ci2ptx: Convert confidence interval to PicTeX format.
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memstack.h"
#include "mismatch.h"
#include "alloc2d.h"
#include "maxlike.h"
#include "bye.h"
#define PROGRAM "CI2PTX"  
/** Dimensions of plots & separations between them in inches: **/
#define PLOTWIDTH 1.0
#define PLOTHEIGHT 0.5
#define HSEP 1.0
#define VSEP 0.75

/******* prototypes *******/
void usage(char *message);
char *pictexval(double pval);
char *stringval(double pval);

/* Main routine */
int    main(int argc, char **argv)
{
  int     i, ntau;
  FILE   *fp=stdout, *ifp=stdin;
  char *fname;
  double plotheight = PLOTHEIGHT;
  double plotwidth = PLOTWIDTH;
  /* indices for next line:0=theta0,1=theta1,2=tau */
  double from[3], to[3], by[3], needrange[3]; 
  double log10theta0, growth=0, tau, pval;
  /* for pictex output */
  int   thisfig, figrow, figcol, ntheta0, ntot;
  double  xpoint, ypoint;
  double  minx, maxx, minxtick, maxxtick, xtickvalue[25];
  double  miny, maxy, minytick, maxytick, ytickvalue[25];
  int   nxticks=4, nyticks=4;
  char  **xtickmark, **ytickmark;
  ASSIGNMENT *a;

  fprintf(fp,"%% -*-latex-*-");
/****** command line arguments *********/
  for(i=1; i<argc; i++)
    if(argv[i][0] == '-') switch(argv[i][1])  /* flag args */
    {
    case 'h':
      plotheight = atof(argv[i]+2);
      break;
    case 'w':
      plotwidth = atof(argv[i]+2);
      break;
    case '-':
    case '?':
      usage("");
    default:
      break;
    }else
    {
      fprintf(fp,"\n%%InputFile = %s", argv[i]);
      fflush(fp);
      if(ifp!=stdin)
	usage("Too many input files");
      ifp = (FILE *) mustopen(argv[i], "r");  /* open input file */
      fname = argv[i];
    }
  if(ifp == stdin)
    fname = "stdin";
  
/***** read ranges from data **********/
  needrange[0] = needrange[1] = needrange[2] = 1;
  do{
    a=getassignment(ifp);
    if(a==NULL)
      error("Ran out of data before finding ranges");
    switch(a->lhs)
    {
    case RangeLog10Theta0:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[0] = atof(a->val[0]);
      to[0]   = atof(a->val[1]);
      by[0]   = atof(a->val[2]);
      needrange[0] = 0;
      break;
    case RangeGrowth:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[1] = atof(a->val[0]);
      to[1]   = atof(a->val[1]);
      by[1]   = atof(a->val[2]);
      needrange[1] = 0;
      break;
    case RangeTau:
      if(a->length != 3)
	error("In input data: ranges must have 3 values, from, to, and by");
      from[2] = atof(a->val[0]);
      to[2]   = atof(a->val[1]);
      by[2]   = atof(a->val[2]);
      needrange[2] = 0;
      break;
    default:
      error("Illegal input while reading ranges");
    }
    freeassignment(a);
  }while(needrange[0] || needrange[1] || needrange[2]);

  if(needrange[0])
    error("Need RangeLog10Theta0");
  if(needrange[1])
    error("Need RangeGrowth");
  if(needrange[2])
    error("Need RangeTau");

  fprintf(fp,"\n%%RangeLog10Theta0: %f %f %f ",
	  from[0], to[0], by[0]);
  fprintf(fp,"\n%%RangeGrowth: %f %f %f ",
	  from[1], to[1], by[1]);
  fprintf(fp,"\n%%RangeTau: %f %f %f ",
	  from[2], to[2], by[2]);

  fflush(fp);

  /*** x axis tick marks ***/
  minxtick = minx = from[2];
  maxxtick = maxx = to[2];
  xtickmark = (char**) alloc2d(2*nxticks, 25, sizeof(char));
  if(xtickmark==NULL)
    error("Can't allocate xtickmark");
  lblaxis(&minxtick, &maxxtick, &nxticks, xtickmark, xtickvalue, 0);
  if(minxtick < minx)
    minx = minxtick;        /* min x axis value */
  if(maxxtick >  maxx)
    maxx = maxxtick;        /* max x axis value */

  /*** y axis tick marks ***/
  minytick = miny = from[0];
  maxytick = maxy = to[0];
  ytickmark = (char**) alloc2d(2*nyticks, 25, sizeof(char));
  if(ytickmark==NULL)
    error("Can't allocate ytickmark");
  lblaxis(&minytick, &maxytick, &nyticks, ytickmark, ytickvalue, 1);
  if(minytick < miny)
    miny = minytick;        /* min y axis value */
  if(maxytick > maxy)
    maxy = maxytick;        /* max y axis value */

  fprintf(fp,"\n%%Varying log10theta0 from %g to %g by %g",
	  from[0], to[0], by[0]);
  fprintf(fp,"\n%%Varying log10(theta1/theta0) from %g to %g by %g",
	  from[1], to[1], by[1]);
  fprintf(fp,"\n%%Varying tau from %g to %g by %g",
	  from[2], to[2], by[2]);

  ntau = (to[2] - from[2])/by[2];
  ntheta0 = (to[0] - from[0])/by[0];
  
  /***** Confidence interval *********/
  fprintf(fp,"\n%%%% Confidence Interval: PicTeX output %%%%");
  fprintf(fp,"\n%%%% Definitions:");
  fprintf(fp,"\n\\def\\hunit{%fin}",plotwidth/(maxx - minx));
  fprintf(fp,"\n\\def\\vunit{%fin}",plotheight/(maxy - miny));
  fprintf(fp,"\n\\def\\accept{{$\\bullet$}}");
  fprintf(fp,"\n\\def\\99pct{{$\\cdot$}}");
  fprintf(fp,"\n\\def\\reject{{\\footnotesize$\\circ$}}");
  fprintf(fp,"\n%%%% Use the following to hide the 1%% region");
  fprintf(fp,"\n\\def\\99pct{\\reject}");
  fprintf(fp,"\n%%%%Plots are %f inch wide and %f inches high",
	  plotwidth, plotheight);
  fprintf(fp,"\n\\begin{figure}");
  fprintf(fp,"\n\\begin{center}\n\\framebox{\\mbox{%%\n\\beginpicture");
  fprintf(fp,"\n\\headingtoplotskip=0.5\\baselineskip");
  fprintf(fp,"\n\\valuestolabelleading=0.4\\baselineskip");

  thisfig = figrow = figcol = ntot = 0;
  while( ( a = getassignment(ifp) ) != NULL)
  {
    if(a==NULL)
      break;
    switch(a->lhs)
    {
    case Growth:  /* new value of growth: start new plot */
      growth = atof(a->val[0]);
      if(figcol==0)
	xpoint = maxx + HSEP*(maxx - minx)/plotwidth;
      else
	xpoint = minx;
      ypoint = figrow * (maxy-miny) * (1.0 + VSEP/plotheight) - miny;
      fputs("\n%%%%%%%%%%%% ",fp);
      fprintf(fp, "Plot figure in row %d col %d", figrow, figcol);
      fputs(" %%%%%%%%%%%%",fp);
      fprintf(fp,"\n\\setcoordinatesystem units <%s,%s> point at %f %f",
	      "\\hunit", "\\vunit", xpoint, ypoint);
      fprintf(fp,"\n\\setplotarea x from %f to %f, y from %f to %f",
	      minx, maxx, miny, maxy);
      /** Bottom axis **/
      if(growth == 0)
	fprintf(fp,"\n\\axis bottom shiftedto y=%f /",
		miny - 0.05*(maxy - miny)/plotheight);
      else
      {
	fprintf(fp,"\n\\axis bottom shiftedto y=%f label {$\\tau$}",
		miny - 0.05*(maxy - miny)/plotheight);
	pictex_ticks(fp, nxticks, xtickmark, xtickvalue);
      }
      /** Left axis **/
      fprintf(fp,"\n\\axis left shiftedto x=%f label {$\\theta_0$}",
	      minx - 0.05*(maxx-minx)/plotwidth);
      pictex_ticks(fp, nyticks, ytickmark, ytickvalue);

      if(growth == 0.0)
	fprintf(fp,"\n\\plotheading{\\small No growth}");
      else if(growth == 1.0)
	fprintf(fp,"\n\\plotheading{\\small 10-fold growth}");
      else
	fprintf(fp,"\n\\plotheading{\\small$10^{%g}$-fold growth}",
		growth);

      fprintf(fp,"\n%%%18s %8s %13s", "", "tau", "log10[theta0]");
      thisfig++;
      if(thisfig % 2 == 0)
	figrow++;
      figcol = (figcol == 0);  /* toggle */
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
      log10theta0 = atof(a->val[1]);
      pval = atof(a->val[2]);
      if(pval < 0.0)
      {
	fprintf(fp,"\n%%     %11s   %8.5g %13.5g    p-val=%g",
		"ERROR",  tau, log10theta0, pval);
	continue;
      }
      record_test(tau, log10theta0, growth, pval);
      if(growth == 0)
	tau = 0.5*(maxx + minx);
      fprintf(fp,"\n\\put %11s at %8.5g %13f   %%p-val=%g",
	      pictexval(pval), tau, log10theta0, pval);
      break;
    default:
      error("Illegal lhs reading confidence interval");
    }
    freeassignment(a);
  }
  fprintf(fp,"\n\\endpicture}}\n\\end{center}");
  fprintf(fp,"\n\\caption{Confidence Region}");
  ptxmaxlike(fp);
  fprintf(fp,"\n\\end{figure}");
  putc('\n', fp);
  if(fp != stdout)
    fclose(fp);
  return 0;
}

char *stringval(double pval)
{
  if(pval <= 0.05)
    return("Reject");
  return("Accept");
}

char *pictexval(double pval)
{
  int ip;

  /* re-express as an integer to avoid rounding errors */
  ip = floor(1000.0*pval + 0.5);
  if(ip < 10)
    return("{\\reject}");  /* reject at 1% level */
  if(ip < 50)
    return("{\\99pct}");  /* w/i 99% confidence region */
  return("{\\accept}");   /* w/i 95% confidence region */
}

char *msg = 
"usage:  \"ci2ptx [options] [inputfile]\"\n\
  where options may include\n\
   -h<number>     Set plot height to <number> inches.\n\
   -w<number>     Set plot width to <number> inches.\n\
By default, the program reads standard input.\n";


void usage(char *message)
{
  fprintf(stderr,"\n%s\n%s\n", message, msg);
  exit(1);
}

