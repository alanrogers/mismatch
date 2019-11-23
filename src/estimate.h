#ifndef ESTIMATE_H
#define ESTIMATE_H
/** defined in estimate.c **/
void estimate(SIMULATION *s, int msize, double *theory_f, int *match,
	      unsigned eflags);
void mom2(double *theta0, double *theta1, double *tau,
	  double *k, int *h, int msize);
double getroughness(int *h, int max);
double get_mse(double *f, int *h, int msize);
double get_mae(double *f, int *h, int dim);
#endif /* ESTIMATE_H */
