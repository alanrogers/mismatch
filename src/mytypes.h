#ifndef MYTYPES_H
#define MYTYPES_H
typedef double DOUBLE ; /* for things that may need extra precision */
typedef long LOGICAL, INTEGER;
enum keyword {Semicolon, Equals, Test, Eof, Badval, InputFile, Sampsize,
	      NSites, Histogram, Cumulants, Estimates, Labels,
	      Seg, RangeLog10Theta0, RangeGrowth, RangeTau, Growth};
/* 16 keywords plus Badval */

typedef struct {
  enum keyword lhs;    /* left hand side of equation */
  char      **val;    /* vector of values */
  int          length;
} ASSIGNMENT;

typedef struct test {
  double tau, log10theta0, growth, pval;
} TEST;

typedef struct {
    int lo;          /* f[lo] is lowest permissible index */
    int hi;          /* f[hi] is highest permissible index */
    int chopLo;      /* tmp value of lo */
    int chopHi;      /* tmp value of hi */
    float *f;
    float *b;        /* base of array */
} FloatVec;

typedef struct {
    int lo;          /* f[lo] is lowest permissible index */
    int hi;          /* f[hi] is highest permissible index */
    int chopLo;      /* tmp value of lo */
    int chopHi;      /* tmp value of hi */
    double *f;
    double *b;          /* base of array */
} DblVec;

typedef struct {
    int lo;          /* f[lo] is lowest permissible index */
    int hi;          /* f[hi] is highest permissible index */
    int chopLo;      /* tmp value of lo */
    int chopHi;      /* tmp value of hi */
    int *f;
    int *b;          /* base of array */
} IntVec;



#endif /* MYTYPES_H */
