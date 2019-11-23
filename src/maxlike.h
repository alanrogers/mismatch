#ifndef MAXLIKE_H
#define MAXLIKE_H
void record_test(double tau, double log10theta0, double growth, double pval);
void getmaxlike(double *tau, double *theta0, double *theta1, double *pval);
void prmaxlike(FILE *fp);
void ptxmaxlike(FILE *fp);
#endif /* MAXLIKE_H */
