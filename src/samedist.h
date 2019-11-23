#ifndef SAMEDIST_H
#define SAMEDIST_H
double dist(double *obs1, double *obs2, int dim);
double samedist(double **s1, int n1, double **s2, int n2, int dim,
	       int nrand,
	       double (*d)(double *obs1, double *obs2, int dim) );
#endif /* SAMEDIST_H */
