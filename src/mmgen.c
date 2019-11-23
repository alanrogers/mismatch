#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "memstack.h"
#include "mismatch.h"
#include "bye.h"
#include "unirand.h"
#include "alloc2d.h"
#include "charplot.h"
#include "header.h"
#define PROGRAM "MMGEN"
#define MSIZE 25                /* default size of mismatch distribution */
void getcumulants(int *h, int max, double *m);
void usage(void);
void mixdist(double *df, int msize, double *mf, int nmut);
int nucleotide(int siteval);
void printdata(NODE * node, MUTATION_MODEL mut_model, FILE * ofp);
int savedata(NODE * node, char **mat, int nsubj, int nsite);

#define NQUANT 7

/*** external variables **/
double qval[NQUANT] = { 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 };
int do_cumulants = 0;
int print_cumulants = 0;
int do_depth = 0;
int print_hist = 0;
int do_theory = 0;
int do_mpd = 0;
int do_quantiles = 0;
int do_seg = 0;
int print_seg = 0;
int print_data = 0;
int reset_mutation = 0;         /* should mutation rates be reset each time? */
int iterations = 1;
int n = 10;                     /* sample size */
int print_plot = 0;
int init_n_sites = 100;         /* # of sites in finite sites models */
int verbose = 0;
double shape = 1.0;             /* shape of gamma distribution */
double sum;                     /* sum of mismatch distribution */
extern int *mismatch, ***intermatch;

void
main(int argc, char **argv)
{
    int i, j, msize = MSIZE, mlen = 0, *h, stacksize, seg;
    int subj, site, locus, splitloci;
    char buff[20];
    FILE *fp;
    NODE *tree;
    POPHIST *history;
    SIMULATION simulation, *sim = &simulation;
    double m[3], meanseg;
    double *f, *mf;
    double **x = NULL, *depth_vec = NULL;   /* holds simulated estimates */
    int maxmutcount = 0;        /* maximum mutations allowed */
    double cummean[3], cumvar[3], sd;
    double dmean, dvar, depth;
    int it, ntot;
    double *xaxis, *yaxis;
    double eps, epsp1, unity, mean_mpd;
    double *mpd_vec = NULL, *seg_vec = NULL;
    extern double sim_mpd;      /* mean pairwise diff of simulation */
    MUTATION_MODEL mut_model = infinite;
    char ***data = NULL;

    check_defines();
    /* calculate machine epsilon */
    eps = unity = 1.0;
    do {
        eps *= 0.5;
        epsp1 = eps + unity;
    } while(epsp1 > unity);

    fp = fopen("pophist.ini", "r");

    if(fp == NULL) {
        fprintf(stderr, "\nWarning: Can't open pophist.ini");
        exit(1);
    }
    history = gethistory(fp);
    fclose(fp);
    fp = stdout;                /* write to standard output */

  /****** Print header ***********/
    header(PROGRAM, "(generate simulated data)", fp);
    fprintf(fp, "\n%c Cmd line:", START_COMMENT);
    for(i = 0; i < argc; i++)
        fprintf(fp, " %s", argv[i]);
    putc('\n', fp);

    sim->tree = NULL;
    sim->h = NULL;
    sim->fname = NULL;

    writehistory(fp, history, START_COMMENT);

  /** process command-line arguments **/
    for(i = 1; i < argc; i++)
        if(argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'c':
                do_cumulants = !do_cumulants;
                break;
            case 'C':
                print_cumulants = !print_cumulants;
                if(print_cumulants)
                    do_cumulants = 1;
                break;
            case 'd':
                do_depth = !do_depth;
                break;
            case 'D':
                do_mpd = !do_mpd;
                break;
            case 'g':
                shape = atof(argv[i] + 2);
                break;
            case 'i':
                iterations = atoi(argv[i] + 2);
                break;
            case 'l':
                msize = atoi(argv[i] + 2);  /* set length of mismatch distribution */
                break;
            case 'n':
                n = atoi(argv[i] + 2);
                break;
            case 's':
                init_n_sites = atoi(argv[i] + 2);
                break;
            case 'S':
                do_seg = !do_seg;
                break;
            case 'm':
                switch (argv[i][2]) {
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
                    if(argv[i][3] != '\0') {
                        maxmutcount = atoi(argv[i] + 3);
                        fprintf(stderr, "\nSetting maxmutcount=%d",
                                maxmutcount);
                    }
                    break;
                default:
                    usage();
                }
                break;
            case 'p':
                switch (argv[i][2]) {
                case 'd':
                    print_data = !print_data;
                    break;
                case 'h':
                    print_hist = !print_hist;
                    break;
                case 'p':
                    print_plot = !print_plot;
                    break;
                case 's':
                    print_seg = !print_seg;
                    break;
                }
                break;
            case 'q':
                do_quantiles = !do_quantiles;
                break;
            case 'r':
                reset_mutation = !reset_mutation;
                break;
            case 't':
                do_theory = !do_theory;
                break;
            case 'v':
                verbose = !verbose;
                break;
            default:
                usage();
            }
        } else
            usage();

    if(!do_seg)
        print_seg = 0;

    switch (mut_model) {
    case finite_gamma:
        fprintf(fp,
                "\n%c Mutation model = finite-gamma. sites=%d shape param=%f",
                START_COMMENT, init_n_sites, shape);
        fprintf(fp, "\n%c Mutation rates will be set %s.",
                START_COMMENT,
                (reset_mutation ? "independently for each tree" :
                 "once at the top"));
        break;
    case finite_flat:
        fprintf(fp, "\n%c Mutation model = finite-flat. sites=%d",
                START_COMMENT, init_n_sites);
        break;
    case infinite:
        fprintf(fp, "\n%c Mutation model = infinite sites", START_COMMENT);
        break;
    case stepwise:
        fprintf(fp, "\nSimple stepwise mutational model");
        break;
    default:
        fprintf(stderr, "\nILLEGAL MUTATION MODEL. (Value = %d)",
                (int) mut_model);
        exit(1);
    }

    sim->smpsiz = n;
    sim->nsubs = history->K;    /* number of subdivisions */
    sim->subsiz = (int *) mustalloc(sim->nsubs * sizeof (int));
    /* make sample size evenly divisible by sim->nsubs */
    n = sim->smpsiz = (sim->smpsiz / sim->nsubs) * sim->nsubs;
    for(i = 0; i < sim->nsubs; i++)
        sim->subsiz[i] = sim->smpsiz / sim->nsubs;  /* sizes of subdivisions */

    fprintf(fp, "\n%c total sampsize=%d subdivisions=%d",
            START_COMMENT, sim->smpsiz, sim->nsubs);
    fprintf(fp, "\n%c subdivision sizes:", START_COMMENT);
    for(i = 0; i < sim->nsubs; i++)
        fprintf(fp, " %d", sim->subsiz[i]);
    fflush(stdout);
    stacksize = (2 * n - 1) * sizeof (NODE);
    if(mut_model != infinite)
        stacksize += (2 * n - 1) * init_n_sites * sizeof (STATE);
    stacksize += 4000;
    if(setmemstack(stacksize) == NULL)
        error("Can't set memstack");

    fprintf(fp, "\n%c Using histograms of size %d", START_COMMENT, msize);
    mf = (double *) mustalloc(msize * sizeof (double));
    f = (double *) mustalloc(msize * sizeof (double));
    sim->h = h = (int *) mustalloc(msize * sizeof (int));
    xaxis = (double *) mustalloc(msize * sizeof (double));
    yaxis = (double *) mustalloc(msize * sizeof (double));

    if(print_data && (mut_model == finite_flat || mut_model == finite_gamma)) {
        data = (char ***) alloc3d(iterations, sim->smpsiz, init_n_sites,
                                  sizeof (char));
        if(data == NULL)
            error("main:alloc3d");
    }

    /* get theoretical mismatch distribution */
    hinit_theory(msize, history);
    if(mut_model == stepwise) {

    /***********************************************************
    before calls to mixdist, mf[i] must equal prob that a random
    pair differs by i mutations.  From this, mixdist calculates
    the probability that the pair differs by |i| steps and puts
    this into f.
    ************************************************************/
        mf = f_hist(history, mf, msize);
        mixdist(f, msize, mf, msize);
    } else
        f = fixsum(f_hist(history, f, msize), msize);
    (void) fixsum(f, msize);

/*** initialize random number generator *********/
    initrand(0);

    /*** Set indices *********/
    setindices();

    /*** Set mutation model at top ***/
    init_mutation(mut_model, init_n_sites, shape, reset_mutation);

    /****** echo parameters ************/
    fprintf(fp, "\n%c iterations=%d, sampsize/subpop=%d total sampsize=%d",
            START_COMMENT, iterations, n / (history->K), n);
    if(do_quantiles) {
        fprintf(fp, "\nNote: Quantiles option turns off all other output.");

       /*
	* Below, "#if A" separates one version from another.  I'm not
	* sure which is correct.
	*/
#if A
        x = (double **) alloc2d(dim_s + msize, iterations, sizeof (double));
#else
        x = (double **) alloc2d(dim_s, iterations, sizeof (double));
#endif
        if(x == NULL)
            error("in main(): Can't allocate x");
        fprintf(stderr, "\nCalling multisim");
        if(verbose)
            putc('\n', stderr);
        ntot = multisim(sim, history, iterations, x, f, msize, mut_model);
        fprintf(fp, "\n%c Did %d iterations in all.", START_COMMENT, ntot);
#if A
        for(i = 0; i < (dim_s + msize); i++)
#else
        for(i = 0; i < dim_s; i++)
#endif
            qsort(x[i], (unsigned) iterations, sizeof (double),
                  (int (*)(const void *, const void *)) compar);
        fprintf(fp, "\n\nQUANTILES:");
        fprintf(fp, "\n%6s:", "");
        for(i = 0; i < dim_s; ++i) {
            if(i == ndx_s[i_theta0])
                fprintf(fp, " %11s", "Log10theta0");
            else
                fprintf(fp, " %11s", lbl_s[i]);
        }
#if A
        for(i = 0; i < msize; ++i) {
            sprintf(buff, "h[%d]", i);
            fprintf(fp, " %11s", buff);
        }
#endif
        for(i = 0; i < NQUANT; ++i) {
            fprintf(fp, "\n%6.4f:", qval[i]);
            for(j = 0; j < dim_s; ++j) {
                if(j == ndx_s[i_theta0]) {
                    sd = quantile(qval[i], x[j], iterations);
                    if(sd == 0.0)
                        fprintf(fp, " %11s", "Log_0");
                    else
                        fprintf(fp, " %11.6f", log10(sd));
                } else
                    fprintf(fp, " %11.6f",
                            quantile(qval[i], x[j], iterations));
            }
#if A
            for(j = 0; j < msize; ++j) {
                fprintf(fp, " %11.6f",
                        quantile(qval[i], x[dim_s + j], iterations));
            }
#endif
        }
        putc('\n', fp);
        exit(0);                /* end here if quantiles are requested */
    }

/***** main loop if not doing quantiles **********/
    if(do_depth)
        depth_vec = (double *) mustalloc(iterations * sizeof (double));
    if(do_mpd)
        mpd_vec = (double *) mustalloc(iterations * sizeof (double));
    if(do_seg)
        seg_vec = (double *) mustalloc(iterations * sizeof (double));
    for(i = 0; i < 3; i++)
        cummean[i] = cumvar[i] = 0.0;
    mean_mpd = dmean = dvar = meanseg = 0.0;
    for(it = 0; it < iterations; it++) {
        tree = iscoales(sim->nsubs, sim->subsiz, history, msize, 1.0);
        if(do_seg) {
            meanseg += seg = getsegregating();
            seg_vec[it] = seg;
            if(print_seg)
                fprintf(fp, "\n%d segregating sites in iteration %d", seg,
                        it);
        }

        if(print_hist || print_plot || do_cumulants || do_mpd)
            mlen = getmatch(tree);
        if(do_mpd) {
            mean_mpd += sim_mpd;
            fprintf(fp, "\nmean pairwise diff=%f", sim_mpd);
            mpd_vec[it] = sim_mpd;
        }
        if(print_hist) {
            fprintf(fp, "\nh=");
            sum = 0.0;
            for(i = 0; i < mlen; i++) {
                sum += mismatch[i];
                if((i + 1) % 19 == 0)
                    fprintf(fp, "\n  ");
                fprintf(fp, " %3d", mismatch[i]);
            }
            if(sum > 0.0) {
                fprintf(fp,
                        "\n%c normalized mismatch distribution inPicTeX format:",
                        START_COMMENT);
                fprintf(fp, "\n\\multiput {$\\circ$} at ");
                for(i = 0; i < mlen; i++) {
                    if((i + 1) % 5 == 0)
                        fprintf(fp, "\n  ");
                    fprintf(fp, " %d %f", i, mismatch[i] / sum);
                }
                fprintf(fp, "\n/");
            }
        }
        if(print_plot)
            plotivec(fp, msize, mismatch);
        if(do_cumulants) {
            getcumulants(mismatch, mlen, m);
            for(i = 0; i < 3; i++) {
                cummean[i] += m[i];
                cumvar[i] += m[i] * m[i];
            }
            if(print_cumulants)
                fprintf(fp, "\nm=%f v=%f E[(x-m)^3]=%f", m[0], m[1], m[2]);
        }
        if(do_depth) {
            depth = treedepth(tree);
            depth_vec[it] = depth;
            dmean += depth;
            dvar += depth * depth;
        }

        if(print_data) {
            if(mut_model == finite_flat || mut_model == finite_gamma)
                savedata(tree, data[it], sim->smpsiz, init_n_sites);
            else {
                fputs("\n%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%", fp);
                fprintf(fp, "\n%d %d", n, init_n_sites);
                printdata(tree, mut_model, fp);
                fputs("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", fp);
            }
        }
        stackfree();
    }

    if(do_cumulants) {
        for(i = 0; i < 3; i++) {
            cummean[i] /= iterations;
            if(iterations > 1)
                cumvar[i] =
                    (cumvar[i] -
                     iterations * cummean[i] * cummean[i]) / (iterations -
                                                              1.0);
            else
                cumvar[i] = -9.999;
        }
        fprintf(fp, "\ncummean:");
        for(i = 0; i < 3; i++)
            fprintf(fp, " %15.6g", cummean[i]);
        fprintf(fp, "\n     sd:");
        for(i = 0; i < 3; i++)
            fprintf(fp, " %15.6g", sqrt(cumvar[i]));
    }
    if(do_seg && iterations > 0) {
        meanseg /= iterations;
        fprintf(fp, "\nMean segregating sites: %f", meanseg);
        qsort(seg_vec, (unsigned) iterations, sizeof (double),
              (int (*)(const void *, const void *)) compar);
        fprintf(fp, "\nQUANTILES of segregating sites:");
        for(i = 0; i < NQUANT; i++)
            fprintf(fp, "\n%6.4f: %11.6f",
                    qval[i], quantile(qval[i], seg_vec, iterations));
        putc('\n', fp);
    }
    if(do_mpd && iterations > 0) {
        mean_mpd /= iterations;
        fprintf(fp, "\nMean of mean pairwise diff: %f", mean_mpd);
        qsort(mpd_vec, (unsigned) iterations, sizeof (double),
              (int (*)(const void *, const void *)) compar);
        fprintf(fp, "\nQUANTILES of mpd:");
        for(i = 0; i < NQUANT; i++)
            fprintf(fp, "\n%6.4f: %11.6f",
                    qval[i], quantile(qval[i], mpd_vec, iterations));
        putc('\n', fp);
    }

    if(do_depth && iterations > 1) {
        dmean /= iterations;
        dvar = (dvar - iterations * dmean * dmean) / (iterations - 1.0);
        fprintf(fp, "\ndepthmean: %15.6g", dmean);
        fprintf(fp, "\n depthvar: %15.6g", dvar);
        fprintf(fp, "\n depthSD: %15.6g", sqrt(dvar));
        fprintf(fp, "\n depthSE: %15.6g", sqrt(dvar / iterations));
        qsort(depth_vec, (unsigned) iterations, sizeof (double),
              (int (*)(const void *, const void *)) compar);
        fprintf(fp, "\nQUANTILES of depth:");
        for(i = 0; i < NQUANT; i++)
            fprintf(fp, "\n%6.4f: %11.6f",
                    qval[i], quantile(qval[i], depth_vec, iterations));
        putc('\n', fp);
    }
    if(do_theory > 0) {
        fputs("\nWARNING:", fp);
        fputs(" theoretical mismatch distribution ignores pop structure", fp);
        if(mut_model == stepwise) {
            pictex_rvec(fp, msize, mf, "Pr[random pair differs by i sites]");
            if(print_plot)
                plotrvec(fp, msize, mf);
            pictex_rvec(fp, msize, f, "Pr[random pair differs by |i| steps]");
            if(print_plot)
                plotrvec(fp, msize, f);
        } else {
            pictex_rvec(fp, msize, f, "Pr[random pair differs by i sites]");
            if(print_plot)
                plotrvec(fp, msize, f);
        }
    }

    /* print data matrix in phylip format */
    if(print_data && (mut_model == finite_flat || mut_model == finite_gamma)) {
#ifndef NDEBUG
        assert(data != NULL);
#endif
        splitloci = (init_n_sites > 20 ? 1 : 0);
        fprintf(fp, "\n%d %d", sim->smpsiz, iterations);
        for(locus = 0; locus < iterations; locus++)
            fprintf(fp, " %d", init_n_sites);
        fprintf(fp, "\n%%%9s ", "");
        for(site = 0; site < init_n_sites && site < 50; site++) {
            if(site > 0 && site % 10 == 0)
                putc(' ', fp);
            fprintf(fp, "%1d", site % 10);
        }
        for(subj = 0; subj < sim->smpsiz; subj++) {
            sprintf(buff, "subj%d", subj);
            fprintf(fp, "\n%-10s ", buff);
            for(locus = 0; locus < iterations; locus++) {
                if(locus > 0) {
                    if(splitloci)
                        fprintf(fp, "\n%10s ", "");
                    else {
                        if(init_n_sites > 10) {
                            if(locus * (init_n_sites + 1) == 50)
                                fprintf(fp, "\n%10s ", "");
                            else if(locus > 0 && locus % 10 == 0)
                                putc(' ', fp);
                        }
                    }
                }
                for(site = 0; site < init_n_sites; site++) {
                    if(splitloci) {
                        if(site > 0 && site % 50 == 0)
                            fprintf(fp, "\n%10s ", "");
                        else if(site > 0 && site % 10 == 0)
                            putc(' ', fp);
                    } else {
                        if((locus > 0 || site > 0)
                           && ((locus * init_n_sites + site) % 50 == 0))
                            fprintf(fp, "\n%10s ", "");
                        else if(site > 0 && site % 10 == 0)
                            putc(' ', fp);
                    }
                    fprintf(fp, "%c", nucleotide(data[locus][subj][site]));
                }
            }
        }
    }
    putc('\n', fp);
    exit(0);
}

void
usage(void)
{
    char buff[30];

    fprintf(stderr, "\nusage: mmgen [options]\n  where options may include:");

    option("-c", "Generate cumulants?", YES(do_cumulants));
    option("-C", "Print cumulants?", YES(print_cumulants));
    option("-d", "Calc tree depth?", YES(do_depth));
    sprintf(buff, "%s", YES(do_mpd));
    option("-D", "Calc mean pairwise diff?", YES(do_mpd));
    sprintf(buff, "%f", shape);
    option("-g<x>", "Set gamma shape parameter to <x>", buff);
    sprintf(buff, "%d", iterations);
    option("-i<x>", "Set iterations to <x>.", buff);
    sprintf(buff, "%d", MSIZE);
    option("-l<x>", "Set length of mismatch distribution to <x>.", buff);
    option("-mi", "mutation model = infinite sites", NULL);
    option("-mf", "mutation model = finite sites w/ equal rates", NULL);
    option("-mg", "mutation model = finite sites w/ gamma rates", NULL);
    option("-ms", "mutation model = stepwise", NULL);
    sprintf(buff, "%d", n);
    option("-n<x>", "Set sample size to <x>.", buff);
    option("-pd", "Print data?", YES(print_data));
    option("-ph", "Print histograms?", YES(print_hist));
    option("-pp", "Plot results?", YES(print_plot));
    option("-ps", "Print segregating sites?", YES(print_seg));
    option("-q", "Calc quantiles?", YES(do_quantiles));
    option("-r", "Reset gamma-model mutation rates each time?",
           YES(reset_mutation));
    sprintf(buff, "%d", init_n_sites);
    option("-s<x>", "Set # of sites to <x> (finite sites only)", buff);
    option("-S", "Do segregating sites?", YES(do_seg));
    sprintf(buff, "%s", YES(do_theory));
    option("-t", "Calc theoretical distribution.", buff);
    sprintf(buff, "%s", YES(verbose));
    option("-v", "Verbose mode?", buff);
    putc('\n', stderr);
    exit(1);
}

/****************************************************************
Calculate the probability that two individuals differ by |i| steps.

On entry, mf should be an array of dimension nmut which contains the
theoretical mismatch distribution, i.e. the probability that a pair of
individuals are separated by i mutations.  On return, array df will
contain the probability that the difference between them equals i.
The algorithm is:

         maxmut
  df[0] = sum    mf[y] * prob_step(0,y)
          y=0 

for x=0, and 

         maxmut
  df[x] = sum    mf[y] * (prob_step(x,y) + prob_step(-x,y))
          y=0 

             maxmut
        = 2   sum    mf[y] prob_step(x,y)
              y=0 

if x > 0.
****************************************************************/
void
mixdist(double *df, int msize, double *mf, int nmut)
{
    static double *ff = NULL;
    static ffsize = 0;
    int x, m;

    if(ffsize < msize) {
        if(ff != NULL)
            free(ff);
        ff = (double *) mustalloc(msize * sizeof (double));
        ffsize = msize;
    }
    for(x = 0; x < msize; x++)
        ff[x] = 0.0;
    /* m represents number of mutations */
    for(m = 0; m < nmut; m++) {
        ff[0] += mf[m] * prob_step(0, m);
        /* x represents number of differences */
        for(x = 1; x <= m; x++)
            ff[x] += 2.0 * mf[m] * prob_step(x, m);
    }
    for(x = 0; x < msize; x++)
        df[x] = ff[x];

}

/* translate from internal representation of nucleotides */
int
nucleotide(int siteval)
{
    switch (siteval) {
    case 0:
        return ('A');
    case 1:
        return ('T');
    case 2:
        return ('G');
    case 3:
        return ('C');
    default:
        return ('?');
    }
    return ('?');
}

#if 1

/****************************************************************
savedata: store the data from a tree into a matrix so that it can
all be printed later.  The problem with this routine is that it 
treats the individuals as exchangeable.  If they are not--as for 
example if the population is subdivided--then this routine will 
generate phony data.
****************************************************************/
int
savedata(NODE * node, char **mat, int nsubj, int nsite)
{
    int i, got = 0;

    if(node == NULL)            /* nothing to do */
        return (0);

    if(node->left == NULL && node->right == NULL) { /* leaf */
        if(node->state.seq == NULL)
            error("leaf has no data");
        for(i = 0; i < nsite; i++)
            mat[0][i] = node->state.seq[i];
        return (1);
    }

    /* if we get here, we're in an internal node */
    got += savedata(node->left, mat + got, nsubj - got, nsite);
    got += savedata(node->right, mat + got, nsubj - got, nsite);

    return (got);
}
#endif

void
printdata(NODE * node, MUTATION_MODEL mut_model, FILE * ofp)
{
    int i;
    extern int n_sites;

    if(node == NULL || mut_model == infinite)   /* nothing to print */
        return;
    if(node->left != NULL && node->right != NULL) {
        printdata(node->left, mut_model, ofp);
        printdata(node->right, mut_model, ofp);
        return;
    }
    switch (mut_model) {
    case finite_flat:
    case finite_gamma:
        if(node->state.seq == NULL)
            break;
        fprintf(ofp, "\n%-9x ", (unsigned) node);
        for(i = 0; i < n_sites; i++)
            putc(nucleotide(node->state.seq[i]), ofp);
        break;
    case stepwise:
        fprintf(ofp, "\n%-9x %d", (unsigned) node, node->state.steps);
        break;
    default:
        error("printdata: bad mut_model");
    }
    return;
}
