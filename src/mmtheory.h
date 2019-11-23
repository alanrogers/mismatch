#ifndef MMTHEORY_H
#define MMTHEORY_H
double *getpoisson(double *poi, unsigned dim, double tau);
double *geteq(double *eq, unsigned dim, double theta);
double *nextf(double *f, double *f0, unsigned dim, double theta, double tau);
double *getf(double *f, unsigned dim, double theta0, double theta1,
	     double tau);
double *f_hist(double *f, unsigned dim, POPHIST *ph);
double *f_2param(double *f, unsigned dim, double theta, double tau);
double *fixsum(double *f, unsigned dim);
#endif /* MMTHEORY_H */








