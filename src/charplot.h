#ifndef CHARPLOT_H
#define CHARPLOT_H
/*** Prototypes ***/
void charplot(FILE *fp, int n, double *x, double *y);
int  adjust(double *x, int n, int min, int max, double *xtic, int *nxtic);
int  gettic(char *s, double *tic);
FILE *mustopen(char *name, char *mode);
double round(double x);
void plotrvec(FILE *fp, int n, double *y);
void plotivec(FILE *fp, int n, int *y);
#endif /* CHARPLOT_H */
