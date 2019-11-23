/**************************************************************** 
MKSPEC: reads sequence data in phylip format and writes a site
frequency spectrum.
****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "mismatch.h"
#include "bye.h"
#include "getcic.h"
#include "header.h"
#define NAMELENGTH 10
#define MSIZE 1000
#define BLEN 200
#define MAXPOPS 100

struct data {
    int     subjects;		/* # number of subjects to compare */
    char  **name;
    int     sites;		/* # of nucleotide sites */
    char  **data;		/* data[subject][site] */
};

/****** Prototypes ********/

int     mygetc(FILE * fp);
void  **alloc2d(int dim1, int dim2, int size);
int     free2d(void **mat);
struct data *getdata(FILE * infile);
int     getsubject(char *name, char *data, int sites, FILE * infile);
char   *getname(char *name, FILE * infile);
void    usage(void);
void    makeSpectrum(struct data *dat);

int     curcol;			/* current column */
int     verbose = 0;

int     mygetc(FILE * fp)
{
    int     c;

    c = getcic(fp);		/* getcic ignores comments */
    if (c == '\n')
	curcol = -1;
    else
	++curcol;
    return (c);
}

char   *getname(char *name, FILE * infile)
{
    int     i, c;

    if (infile == NULL)
	return (NULL);
    while (isspace(c = mygetc(infile)) && c != EOF) /* go to start of name */
	;
    if (c == EOF)
	return (NULL);
    if (curcol != 0) {		/* does name start in col 0? */
	fprintf(stderr, "\nData are incorrectly formatted.");
	fprintf(stderr, "  Name appears to begin in column=%d w/ char %c.",
		curcol, c);
	fprintf(stderr, "\n  Should have been in col 0.\n");
	exit(1);
    }
    i = 0;
    name[i++] = c;
    while (curcol < 9) {
	c = mygetc(infile);
	if (c == EOF) {
	    name[i] = '\0';
	    return (name);
	}
	if (isspace(c))
	    c = '_';		/* blanks become underscores */
	name[i++] = c;
    }
    while (name[i - 1] == '_')
	i--;			/* strip trailing underscores */
    name[i] = '\0';
    if (i > 0)
	return (name);		/* normal return */
    return (NULL);		/* error: return NULL if name is empty */
}

struct data *getdata(FILE * infile)
{
    int     i, j, c, first_subject;
    char    buff[BLEN];
    char   *ref;
    const int use_ref = 1;
    struct data *s;

    if (infile == NULL)
	error("getdata was handed a null file pointer");

    s = (struct data *) mustalloc(sizeof(struct data));

    /* get numbers of subjects and sites */
    if (getintic(&(s->subjects), buff, BLEN, infile) == EOF)
	return (NULL);
    if (getintic(&(s->sites), buff, BLEN, infile) == EOF)
	return (NULL);
    fflush(stdout);
    fprintf(stderr,
	    "\n%c 1st line of input specifies %d subjects and %d sites",
	    START_COMMENT, s->subjects, s->sites);
    /* skip rest of line */
    do {
	c = getcic(infile);
    } while (c != '\n' && c != EOF);

    curcol = -1;
    if (verbose)
	printf("\n%c Expecting %d subjects X %d sites = %d sites in all\n",
	       START_COMMENT,
	       s->subjects, s->sites, (s->subjects) * (s->sites));
    if (s->subjects == 0 || s->sites == 0)
	return (NULL);
    s->name = (char **) alloc2d(s->subjects, NAMELENGTH + 1, sizeof(char));
    if (s->name == NULL)
	error("out of memory in getdata()");
    s->data = (char **) alloc2d(s->subjects, s->sites, sizeof(char));
    if (s->data == NULL)
	error("out of memory in getdata()");
    ref = (char *) mustalloc(s->sites * sizeof(char));

  /** read reference line **/
    if (getsubject(s->name[0], ref, s->sites, infile) < s->sites)
	error("getdata: EOF reading reference line");
    for (j = 0; j < s->sites; j++)
	if (ref[j] == '.') {
	    fprintf(stderr, "\nName on reference line: \"%s\"",
		    s->name[0]);
	    fprintf(stderr, "\nData on reference line: \"%s\"",
		    s->name[0]);
	    for (i = 0; i < s->sites; i++)
		fprintf(stderr, "%c", ref[i]);
	    error
		("The \".\" character is not allowed in the reference line");
	}

    if (use_ref) {		/* reference line is part of data */
	for (j = 0; j < s->sites; j++)
	    s->data[0][j] = ref[j];
	first_subject = 1;
    } else			/* reference line is not part of data */
	first_subject = 0;

    for (i = first_subject; i < s->subjects; i++) {
	if (getsubject(s->name[i], s->data[i], s->sites, infile) <
	    s->sites) {
	    fprintf(stderr, "\ngetdata: Insufficient data for subject %d.",
		    i);
	    exit(1);
	}
	/* replace '.' w/ corresponding value from reference sequence */
	for (j = 0; j < s->sites; j++)
	    if (s->data[i][j] == '.')
		s->data[i][j] = ref[j];
    }
    i = 0;
    while ((c = getcic(infile)) != EOF)
	if (!isspace(c)) {
	    if (i == 0) {
		fprintf(stderr,
			"\nExtraneous characters at end of data:\n");
		i = 1;
	    }
	    putc(c, stderr);
	}
    if (i == 1)
	exit(1);
    return (s);
}

/** read data for one subject; return number of sites read **/
int     getsubject(char *name, char *data, int sites, FILE * infile)
{
    int     j, c;

    if (getname(name, infile) == NULL)
	error("getsubject: subject has no name.");
    for (j = 0; j < sites; j++) {
	while (isspace(c = mygetc(infile)))	/* skip blanks */
	    ;
	if (c == EOF)
	    return (j);
	data[j] = toupper(c);
    }
    return (sites);
}

/****************************************************************************
I normally keep the following programs in a separate file called "alloc2d.c",
but have concatenated everything for ease of distribution.
****************************************************************************/
/*
 * Matrix allocation routines.  
 * To allocate a 10 by 10 matrix of doubles use:
 *
 *	void **alloc2d();
 *	double **x;
 *
 *	x = (double **) alloc2d(10, 10, sizeof(double));
 *
 * To free this matrix use:
 *
 *	free2d(x);
 */
/*--Allocate & free 2 arrays--*/
void  **alloc2d(int dim1, int dim2, int size)
{
    int     i;
    unsigned nelem;
    char   *p, **pp;

    nelem = dim1 * dim2;
    p = (void *) malloc((unsigned) nelem * size);
    if (p == NULL)
	return (NULL);
    pp = (char **) malloc((unsigned) dim1 * (unsigned) sizeof(char *));
    if (pp == NULL) {
	free(p);
	return (NULL);
    }

    for (i = 0; i < dim1; i++)
	pp[i] = p + i * dim2 * size;

    return ((void **) pp);
}
int     free2d(void **mat)
{
    if (mat != NULL && *mat != NULL)
	free((void *) *mat);
    if (mat != NULL)
	free((void *) mat);
    return (0);
}

void    usage(void)
{
    fflush(stdout);
    fprintf(stderr, "\nUsage: mkspec [filename1 ...]");
    fprintf(stderr, "\n\nFile format:");
    fprintf(stderr, "\n#_of_subjects #_of_sites");
    fprintf(stderr, "\nReference ATGATGATG");
    fprintf(stderr, "\nName1_____.....A.A.");
    fprintf(stderr, "\nName2_____G.C......");
    fprintf(stderr, "\nName3_____.A.GA.G.C\n");
    fprintf(stderr, "\nEtc.");
    fprintf(stderr, "\nDots indicate identity w/ reference sequence.\n");
    fprintf(stderr,
	    "\nMultiple input files generate intermatch comparisons\n");
    exit(1);
}

/* make mismatch or intermatch distributions */
void    makeSpectrum(struct data *dat)
{
    int    *spectrum;
    int     segregating=0;  /* number of segregating sites */
    int     i, j;
    int     nA, nT, nG, nC, nMinor, nucleotide;
    int     max;
    double  a, e, theta;

    /* max is maximal count of minor allele */
    if (dat->subjects % 2 == 0)	/* even number of sequences */
	max = dat->subjects / 2;
    else			/* odd number of sequences */
	max = (dat->subjects - 1) / 2;

    /*
     * Legal values of spectrum are spectrum[0]..spectrum[max],
     * but we will ignore the 0'th entry.
     */
    spectrum = (int *) malloc((max + 1) * sizeof(int));

    for (i = 0; i < dat->sites; ++i) {
	nA = nT = nG = nC = 0;
	for (j = 0; j < dat->subjects; ++j) {
	    if (dat->data[j][i] == '.')
		nucleotide = dat->data[0][i];
	    else
		nucleotide = dat->data[j][i];
	    switch (toupper(nucleotide)) {
	    case 'A':
		++nA;
		break;
	    case 'T':
		++nT;
		break;
	    case 'G':
		++nG;
		break;
	    case 'C':
		++nC;
		break;
	    default:
		fprintf(stderr, "\nIllegal value for nucleotide: %c\n",
			nucleotide);
		exit(1);
	    }
	}
	/* A */
	nMinor = nA;

	/* T */
	if (nMinor == 0)
	    nMinor = nT;
	else if (nT > 0 && nT < nMinor)
	    nMinor = nT;

	/* G */
	if (nMinor == 0)
	    nMinor = nG;
	else if (nG > 0 && nG < nMinor)
	    nMinor = nG;

	/* C */
	if (nMinor == 0)
	    nMinor = nC;
	else if (nC > 0 && nC < nMinor)
	    nMinor = nC;

	/* Now nMinor = count of least frequent nucleotide, */
	/* excluding zeroes */

	if (nMinor == dat->subjects)	/* ignore monomorphic sites */
	    continue;

	++segregating;        /* count segregating sites */
	++spectrum[nMinor];   /* accumulate spectrum */
    }

    a = 0.0;
    for(i=1; i < dat->subjects; ++i)
	a += 1.0/i;

    theta = segregating/a;

    a = 0.0;
    putchar('\n');
    printf("%6s %6s %8s\n", "Count", "", "Expected");
    printf("%6s %6s %8s\n", "  of", "Number", "number");
    printf("%6s %6s %8s\n", "minor", "of", "of");
    printf("%6s %6s %8s\n", "allele", "sites", "sites");
    printf("%6s %6s %8s\n", "------", "------", "------");
    for (i = 1; i <= max; ++i) {
	if(i == dat->subjects -i)
	    e = theta/i;
        else   /* for folded spectrum */
	    e = theta*(1.0/i + 1.0/(dat->subjects - i));
	a += e;
	printf("%6d %6d %8g\n", i, spectrum[i], e);
    }
    putchar('\n');
    free(spectrum);
}

int     main(int argc, char **argv)
{
    int     i;
    struct data *dat = NULL;
    char   *fname;
    FILE   *fp;

    header("mkspec", "(make a site frequency spectrum)", stdout);
    echo_cmdline(argc, argv);

    fp = stdin;
    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {
	    default:
		usage();
	} else {
	    if (dat != NULL)
		usage();
	    fname = argv[i];
	    fp = (FILE *) fopen(fname, "r");
	    if (fp == NULL) {
		fprintf(stderr,
			"\nError:Can't open input file \"%s\".\n", fname);
		exit(1);
	    }
	    dat = getdata(fp);
	}
    }

    if (dat == NULL)
	error("No data");

    makeSpectrum(dat);

    fflush(stdout);
    putc('\n', stderr);
    return 0;
}
