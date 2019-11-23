#ifndef EPRINTF_H
#define EPRINTF_H
void eprintf(char *progname, char *fmt, ...);
char *estrdup(char *progname, char *s);
void *emalloc(char *progname, size_t n);
#endif /* EPRINTF_H */
