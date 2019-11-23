/**
 * seqstat: reads sequence data in phylip format and writes various
 * descriptive statistics.
 */

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
/* next 3 defs refer to the style in which output is produced */
#define MMEST      0
#define NORMALIZED 1
#define PICTEX     2

struct data {
    int         subjects;       /* # number of subjects to compare */
    char      **name;
    int         sites;          /* # of nucleotide sites */
    char      **data;           /* data[subject][site] */
};

/****** Prototypes ********/
int         mygetc(FILE * fp);
void      **alloc2d(int dim1, int dim2, int size);
int         free2d(void **mat);
struct data *getdata(FILE * infile);
int         getsubject(char *name, char *data, int sites, FILE * infile);
double      makemm(struct data *s1, struct data *s2, char *hit);
int         makeFoldedSpectrum(struct data *s, IntVec * spec,
                               int *minorCount);
char       *getname(char *name, FILE * infile);
void        usage(void);
void        printPolymorphicSites(struct data *s, char *hit);

/** global variables **/
int         curcol;             /* current column */
int         use_ref = 1;        /*is reference line part of data? */
int         out_style = MMEST;
int         printMinor = 0;     /* print counts of minor alleles? */
int         verbose = 0;
int         echoSeq = 0;        /* print the sequence data? */

int mygetc(FILE * fp) {
    int         c;

    c = getcic(fp);             /* getcic ignores comments */
    if(c == '\n')
        curcol = -1;
    else
        ++curcol;
    return (c);
}

char       *getname(char *name, FILE * infile) {
    int         i, c;

    if(infile == NULL)
        return (NULL);

    /* go to start of name */
    while(isspace(c = mygetc(infile)) && c != EOF) ;
    if(c == EOF)
        return (NULL);
    if(curcol != 0) {           /* does name start in col 0? */
        fprintf(stderr, "\nData are incorrectly formatted.");
        fprintf(stderr, "  Name appears to begin in column=%d w/ char %c.",
                curcol, c);
        fprintf(stderr, "\n  Should have been in col 0.\n");
        exit(1);
    }
    i = 0;
    name[i++] = c;
    while(curcol < 9) {
        c = mygetc(infile);
        if(c == EOF) {
            name[i] = '\0';
            return (name);
        }
        if(isspace(c))
            c = '_';            /* blanks become underscores */
        name[i++] = c;
    }
    while(name[i - 1] == '_')
        i--;                    /* strip trailing underscores */
    name[i] = '\0';
    if(i > 0)
        return (name);          /* normal return */
    return (NULL);              /* error: return NULL if name is empty */
}

struct data *getdata(FILE * infile) {
    int         i, j, c, first_subject;
    char        buff[BLEN];
    char       *ref;
    struct data *s;

    if(infile == NULL)
        error("getdata was handed a null file pointer");

    s = (struct data *) mustalloc(sizeof(struct data));

    /* get numbers of subjects and sites */
    if(getintic(&(s->subjects), buff, BLEN, infile) == EOF)
        return (NULL);
    if(getintic(&(s->sites), buff, BLEN, infile) == EOF)
        return (NULL);
    fflush(stdout);
    fprintf(stderr,
            "\n%c 1st line of input specifies %d subjects and %d sites",
            START_COMMENT, s->subjects, s->sites);
    /* skip rest of line */
    do {
        c = getcic(infile);
    } while(c != '\n' && c != EOF);

    curcol = -1;
    if(verbose)
        printf("\n%c Expecting %d subjects X %d sites = %d sites in all\n",
               START_COMMENT,
               s->subjects, s->sites, (s->subjects) * (s->sites));
    if(s->subjects == 0 || s->sites == 0)
        return (NULL);
    s->name = (char **) alloc2d(s->subjects, NAMELENGTH + 1, sizeof(char));
    if(s->name == NULL)
        error("out of memory in getdata()");
    s->data = (char **) alloc2d(s->subjects, s->sites, sizeof(char));
    if(s->data == NULL)
        error("out of memory in getdata()");
    ref = (char *) mustalloc(s->sites * sizeof(char));

/** read reference line **/
    if(getsubject(s->name[0], ref, s->sites, infile) < s->sites)
        error("getdata: EOF reading reference line");
    for(j = 0; j < s->sites; j++)
        if(ref[j] == '.') {
            fprintf(stderr, "\nName on reference line: \"%s\"", s->name[0]);
            fprintf(stderr, "\nData on reference line: \"%s\"", s->name[0]);
            for(i = 0; i < s->sites; i++)
                fprintf(stderr, "%c", ref[i]);
            error("The \".\" character is not allowed in the reference line");
        }
    if(use_ref) {               /* reference line is part of data */
        for(j = 0; j < s->sites; j++)
            s->data[0][j] = ref[j];
        first_subject = 1;
    } else                      /* reference line is not part of data */
        first_subject = 0;

    for(i = first_subject; i < s->subjects; i++) {
        if(getsubject(s->name[i], s->data[i], s->sites, infile) < s->sites) {
            fprintf(stderr, "\ngetdata: Insufficient data for subject %d.",
                    i);
            exit(1);
        }
        /* replace '.' w/ corresponding value from reference sequence */
        for(j = 0; j < s->sites; j++)
            if(s->data[i][j] == '.')
                s->data[i][j] = ref[j];
    }
    i = 0;
    while((c = getcic(infile)) != EOF)
        if(!isspace(c)) {
            if(i == 0) {
                fprintf(stderr, "\nExtraneous characters at end of data:\n");
                i = 1;
            }
            putc(c, stderr);
        }
    if(i == 1)
        exit(1);
    return (s);
}

/** read data for one subject; return number of sites read **/
int getsubject(char *name, char *data, int sites, FILE * infile) {
    int         j, c;

    if(getname(name, infile) == NULL)
        error("getsubject: subject has no name.");
    for(j = 0; j < sites; j++) {
        while(isspace(c = mygetc(infile)))  /* skip blanks */
            ;
        if(c == EOF)
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
 *      void **alloc2d();
 *      double **x;
 *
 *      x = (double **) alloc2d(10, 10, sizeof(double));
 *
 * To free this matrix use:
 *
 *      free2d(x);
 */
/*--Allocate & free 2 arrays--*/
void      **alloc2d(int dim1, int dim2, int size) {
    int         i;
    unsigned    nelem;
    char       *p, **pp;

    nelem = dim1 * dim2;
    p = (void *) malloc((unsigned) nelem * size);
    if(p == NULL)
        return (NULL);
    pp = (char **) malloc((unsigned) dim1 * (unsigned) sizeof(char *));
    if(pp == NULL) {
        free(p);
        return (NULL);
    }
    for(i = 0; i < dim1; i++)
        pp[i] = p + i * dim2 * size;

    return ((void **) pp);
}
int free2d(void **mat) {
    if(mat != NULL && *mat != NULL)
        free((void *) *mat);
    if(mat != NULL)
        free((void *) mat);
    return (0);
}

void usage(void) {
    fflush(stdout);
    fprintf(stderr,
            "\nUsage: seqstat [options] [filename1 ...]\n where options may include:");
    option("-c", "print counts of minor alleles", YES(printMinor));
    option("-e", "echo sequence data?", YES(echoSeq));
    option("-r", "is reference line a subject too?", YES(use_ref));
    option("-m", "mmest style output.", YES(out_style == MMEST));
    option("-n", "normalized histograms", YES(out_style == NORMALIZED));
    option("-p", "PicTeX-style output", YES(out_style == PICTEX));
    option("-v", "verbose output", "No");
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
double makemm(struct data *s1, struct data *s2, char *hit) {
    int       **xtab;           /* cross tabulation */
    int        *h;              /* histogram */
    int         i, j, k, msize = 0;
    double      diff, meandiff, sum;

#ifndef NDEBUG
    assert(s1->sites == s2->sites);
#endif
    xtab = (int **) alloc2d(s1->subjects, s2->subjects, sizeof(int));
    if(xtab == NULL)
        error("makemm: can't allocate subjectXsubject matrix");
    for(i = 0; i < s1->subjects; i++)
        for(j = 0; j < s2->subjects; j++) {
            if(s1 == s2 && j >= i)  /* lower triangle only */
                continue;
            xtab[i][j] = 0;
        }

    meandiff = 0.0;
    for(i = 0; i < s1->subjects; i++) {
        for(j = 0; j < s2->subjects; j++) {
            if(s1 == s2 && j >= i)  /* lower triangle only */
                continue;
            diff = 0;
            for(k = 0; k < s1->sites; k++)
                if(s1->data[i][k] != s2->data[j][k]) {
                    hit[k] = 1;
                    diff++;
                }
            xtab[i][j] = diff;
            if(msize <= diff)
                msize = diff + 1;
            meandiff += diff;
        }
    }
    if(s1 == s2)
        meandiff /= s1->subjects * (s2->subjects - 1) / 2;
    else
        meandiff /= s1->subjects * s2->subjects;

    /* allocate array for histogram */
    h = (int *) malloc(msize * sizeof(int));
    if(h == NULL) {
        fflush(stdout);
        fprintf(stderr, "\nmakemm: no memory\n");
        exit(1);
    }
    /* calculate histogram */
    for(i = 0; i < msize; i++)
        h[i] = 0;
    for(i = 0; i < s1->subjects; i++)
        for(j = 0; j < s2->subjects; j++) {
            if(s1 == s2 && j >= i)  /* lower triangle only */
                continue;
            if(xtab[i][j] < msize)
                h[xtab[i][j]] += 1;
            else
                h[msize - 1] += 1;
        }
    while(msize > 0 && h[msize - 1] == 0)
        --msize;

    /* calculate sum */
    sum = 0.0;
    for(i = 0; i < msize; i++)
        sum += h[i];

    /* produce output */
    switch (out_style) {
    case MMEST:                /* output in style needed by mmest */
        printf("\n    mismatch =");
        for(i = 0; i < msize; i++)
            printf(" %d", h[i]);
        printf(" ;");
        break;
    case NORMALIZED:           /* write normalized histogram */
        printf("\n\n    %cNORMALIZED \n    mismatch =", START_COMMENT);
        for(i = 0; i < msize; i++) {
            if((i + 1) % 5 == 0)
                fputs("\n    ", stdout);
            printf(" %f", (double) (h[i] / sum));
        }
        fputs(" ;", stdout);
        break;
    case PICTEX:               /* PicTeX format */
        printf("\n\n    %cPicTeX FORMAT:", START_COMMENT);
        printf("\n    \\plot\n    ");
        for(i = 0; i < msize; i++) {
            if((i + 1) % 5 == 0)
                fputs("\n    ", stdout);
            printf(" %d %f", i, (double) h[i] / sum);
        }
        printf("\n    /\n");
        break;
    default:
        error("Bad value in switch");
    }
    free(h);
    free2d((void **) xtab);
    return (meandiff);
}

/*
 * Make folded site frequency spectrum.
 * On return, spec[i] = number of sites whose minor allele is present
 * in i+1 copies.  The indices of IntVec spectrum should range from
 * 1 through floor(sites/2).  Function returns the number of
 * segregating sites.
 */
int makeFoldedSpectrum(struct data *s, IntVec * spec, int *minorCount) {
    int         site, seq, nA, nT, nG, nC, min;
    int         segregatingSites = 0;

    /* initialize vectors */
    initIntVec(spec);
    memset(minorCount, 0, s->sites * sizeof(int));

    for(site = 0; site < s->sites; site++) {
        nA = nT = nG = nC = 0;
        for(seq = 0; seq < s->subjects; seq++) {
            switch (s->data[seq][site]) {
            case 'a':
            case 'A':
                nA++;
                break;
            case 't':
            case 'T':
                nT++;
                break;
            case 'g':
            case 'G':
                nG++;
                break;
            case 'c':
            case 'C':
                nC++;
                break;
            default:
                fprintf(stderr, "Bad nucleotide: %c in seq %d at site %d\n",
                        s->data[seq][site], seq, site);
                exit(1);
            }
        }
        min = s->subjects;
        if(nA > 0 && nA < min)
            min = nA;
        if(nT > 0 && nT < min)
            min = nT;
        if(nG > 0 && nG < min)
            min = nG;
        if(nC > 0 && nC < min)
            min = nC;
        if(min == s->subjects)
            continue;
        segregatingSites++;
        assert(min <= spec->hi);
        assert(min >= spec->lo);
        spec->f[min]++;
        minorCount[site] = min;
    }
    /* chop off terminal zeroes */
    for(spec->chopHi = spec->hi; spec->f[spec->chopHi] == 0; spec->chopHi--) ;
    return (segregatingSites);
}

/*
 * print polymorphic sites
 */
void printPolymorphicSites(struct data *s, char *hit) {
    int         seq, site, psite;
    assert(s != NULL);
    assert(s->name != NULL);
    assert(s->name[0] != NULL);
    assert(s->data != NULL);
    assert(s->data[0] != NULL);

    /* print reference sequence */
    printf("\n%-10s", s->name[0]);
    for(psite = site = 0; site < s->sites; site++) {
        if(hit[site] == 0)
            continue;
        if(psite > 0 && psite % 50 == 0)
            printf("\n%-10s", "");
        else if(psite > 0 && psite % 10 == 0)
            putchar(' ');
        putchar(s->data[0][site]);
        psite++;
    }

    /* print the other sequences */
    for(seq = 1; seq < s->subjects; seq++) {
        printf("\n%-10s", s->name[seq]);
        for(psite = site = 0; site < s->sites; site++) {
            if(hit[site] == 0)
                continue;
            if(psite > 0 && psite % 50 == 0)
                printf("\n%-10s", "");
            else if(psite > 0 && psite % 10 == 0)
                putchar(' ');
            if(s->data[seq][site] == s->data[0][site])
                putchar('.');
            else
                putchar(s->data[seq][site]);
            psite++;
        }
    }
}

int main(int argc, char **argv) {
    enum { NCOLS = 4 };
    int         row, nrows;
    int         i, j, pop;
    struct data *s[MAXPOPS];
    int         npop = 0;       /* how many populations? */
    int         crease;         /* where to fold spectrum */
    char       *hit;
    int         seg;
    int        *minorCount;     /* counts of minor allele at each site */
    IntVec     *spectrum = NULL;
    char       *fname;
    double      a, meanDiff;

    FILE       *fp;

    header("seqstat", "(descriptive statistics from sequence data)", stdout);
    echo_cmdline(argc, argv);

    fp = stdin;
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-')
            switch (argv[i][1]) {
            case 'c':
                printMinor = !printMinor;
                break;
            case 'e':
                echoSeq = !echoSeq;
                break;
            case 'm':
                out_style = MMEST;
                break;
            case 'n':
                out_style = NORMALIZED;
                break;
            case 'p':
                out_style = PICTEX;
                break;
            case 'r':
                use_ref = TOGGLE(use_ref);
                break;
            case 'v':
                verbose = 1;
                break;
            default:
                usage();
        } else {
            fname = argv[i];
            fp = (FILE *) fopen(fname, "r");
            if(fp == NULL) {
                fprintf(stderr,
                        "\nError:Can't open input file \"%s\".\n", fname);
                exit(1);
            }
            s[npop] = getdata(fp);
            if(s[npop] == NULL)
                error("No data");
            npop += 1;

            /* make sure all pops have same number of sites */
            if(npop > 1 && s[npop - 1]->sites != s[0]->sites) {
                fflush(stdout);
                fprintf(stderr,
                        "\nNumber of sites in pop %d doesn't match pop 0.\n",
                        npop - 1);
                exit(1);
            }
        }
    }

    if(npop == 0)
        error("No data");

    fprintf(stdout,
            "\n%c Results %s include the reference sequence (line 1).",
            START_COMMENT, (use_ref ? "will" : "will not"));

    hit = (char *) mustalloc(s[0]->sites * sizeof(char));
    minorCount = (int *) mustalloc(s[0]->sites * sizeof(int));
    memset(hit, 0, (size_t) s[0]->sites * sizeof(char));

    for(pop = 0; pop < npop; pop++) {
        /* mismatch distributions */
        printf("\n\n%c Population %d", START_COMMENT, pop);
        printf("\n    sequences       : %d", s[pop]->subjects);
        printf("\n    sites           : %d", s[pop]->sites);

        meanDiff = makemm(s[pop], s[pop], hit);
        printf("\n    meanPairwiseDiff");
        printf("\n        per sequence: %g", meanDiff);
        printf("\n        per site    : %g", meanDiff / s[pop]->sites);

        /* spectra */
        crease = (s[pop]->subjects) / 2;
        if(spectrum != NULL && spectrum->hi != crease) {
            free(spectrum->b);
            free(spectrum);
        }
        spectrum = newIntVec(1, crease);

        seg = makeFoldedSpectrum(s[pop], spectrum, minorCount);
        printf("\n   segregating sites: %d", seg);

        /* sum goes from small numbers to big numbers for */
        /* for greater accuracy */
        a = 0.0;
        for(i = s[pop]->subjects - 1; i >= 1; i--)
            a += 1.0 / i;

        printf("\n   theta estimated from segregating sites");
        printf("\n        per sequence: %g", seg / a);
        printf("\n        per site    : %g", seg / (a * s[pop]->sites));

        switch (out_style) {
        case MMEST:
            printf("\n    spectrum = ");
            for(i = spectrum->lo; i <= spectrum->chopHi; i++) {
                if(i > 1 && (i - 1) % 20 == 0)
                    fputs("\n    ", stdout);
                printf(" %d", spectrum->f[i]);
            }
            fputs(" ;", stdout);
            break;
        case NORMALIZED:       /* normalized spectrum */
            printf("\n\n%c   NORMALIZED \n    spectrum =", START_COMMENT);
            for(i = spectrum->lo; i <= spectrum->chopHi; i++) {
                if(i > 1 && (i - 1) % 5 == 0)
                    fputs("\n    ", stdout);
                printf(" %f", ((double) spectrum->f[i]) / seg);
            }
            fputs(" ;", stdout);
            break;
        case PICTEX:           /* PicTeX format */
            printf("\n\n    %c PicTeX FORMAT:", START_COMMENT);
            printf("\n    \\plot\n ");
            for(i = spectrum->lo; i <= spectrum->chopHi; i++) {
                if(i > 1 && (i - 1) % 20 == 0)
                    fputs("\n    ", stdout);
                printf(" %d %d", i, spectrum->f[i]);
            }
            printf("\n    /\n");
            break;
        default:
            error("main: Bad value in switch");
        }

        /* site-by-site statistics */
        if(printMinor) {
            printf("\n    %c Count of minor allele at each site",
                   START_COMMENT);
            printf(" (count=0 for monomorphic sites):");
            nrows = s[pop]->sites / NCOLS;
            if(s[pop]->sites % NCOLS > 0)
                nrows += 1;

            putchar('\n');
            for(i = 0; i < NCOLS; i++) {
                if(i > 0)
                    fputs(" |", stdout);
                printf("%5s %5s", "site", "count");
            }

            for(row = 0; row < nrows; row++) {
                putchar('\n');
                for(i = row; i < s[pop]->sites; i += nrows) {
                    if(i > row)
                        fputs(" |", stdout);
                    printf("%5d %5d", i, minorCount[i]);
                }
            }
        }
        /* echo sequence data */
        if(echoSeq) {
            printf("\n\n%c Polymorphic sites from pop %d:",
                   START_COMMENT, pop);
            printPolymorphicSites(s[pop], hit);
        }
    }

    /* intermatch distributions */
    if(npop > 1)
        printf("\n\n%c Intermatch distributions", START_COMMENT);
    for(i = 1; i < npop; i++)
        for(j = 0; j < i; j++) {
            printf("\n\n%c Populations %d X %d", START_COMMENT, i, j);
            makemm(s[i], s[j], hit);
        }

    fflush(stdout);
    putc('\n', stderr);
    exit(0);
}
