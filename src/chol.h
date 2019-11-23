#ifndef CHOL_H
#define CHOL_H
typedef double real;
int   cholesky(real **A, int n);
void  matmultLLT(real **L, real **B, int n);
real  log_determinant(real **L, int n);
real  dotprod(real *x, real *y, int n);
int   col_L_solve(real **L, int n, real *b);
int   row_L_solve(real **L, int n, real *b);
real  mahal(real *x, real **L, int n);
void  error(char *msg);
void  cholinv(real **L, real **s, int n);
#endif /* CHOL_H */
