#ifndef GETCIC_H
#define GETCIC_H

int   getcic(FILE *fp);
char *getwordic(char *buff, int bufsiz, FILE * ifp);
int   getdoubleic(double *x, char *buff, int bufsiz, FILE *ifp);
int   getintic(int *i, char *buff, int bufsiz, FILE *ifp);

#endif /* GETCIC_H */
