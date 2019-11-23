#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chol.h"
/*****************************************************************
 * These routines compute and manipulate the Cholesky decomposition of
 * a positive definite matrix.
 *   The Cholesky factor is a lower triangular matrix, L, such that
 * L * L' = A, where the "'" denotes matrix transposition.
 *
 * Alan R. Rogers, Dept. of Anthropology, University of Utah, S.L.C.,
 * UT 84108.    February 23, 1998
****************************************************************/
#if 0
/********** Prototypes **********************/
int   cholesky(real **A, int n);
real log_determinant(real **L, int n);
real dotprod(real *x, real *y, int n);
void  matmultLLT(real **L, real **B, int n);
int   col_L_solve(real **L, int n, real *b);
int   row_L_solve(real **L, int n, real *b);
real mahal(real *x, real **L, int n);
#endif

/*****************************************************************
 * Find the left Cholesky factor, L, of a positive definite
 * matrix, A.  The code here was modified from the gaxpy version of
 * the cholesky algorithm on p. 143 of Matrix Computations (2nd ed),
 * by Golub and Van Loan. 
 * 
 * On entry
 *     A     is an n X n matrix of reals as allocated by alloc2d().
 *     n     is the dimension of A.
 * On return
 *     A     contains in its lower triangle the left (lower) Cholesky
 *           factor, L, of A.  L is a lower triangular matrix such
 *           that L L' = A.
 * The routine returns 0 if successful, or 1 if A is not positive 
 * definite.
 *                   Alan R. Rogers Feb 23, 1998
 *****************************************************************/
int cholesky(real **A, int n)
{
  int j,k;
  real root;

  for(j=0; j<n; j++)
  {
    if(j>0)
    {
      for(k=j; k<n; k++)
	A[k][j] -= dotprod(A[k], A[j], j);
    }
#if 1
    if(A[j][j] <= 0.0)
      return(1);
#else
    if(A[j][j] < FLT_EPSILON)
    {
      printf("\nWarning: cholesky fixed up a singular matrix");
      A[j][j] = FLT_EPSILON;
    }
#endif
    root = sqrt(A[j][j]);
    for(k=j; k<n; k++)
      A[k][j] /= root;
  }
  return(0);
}

/* multiply a lower triangular matrix L by its transpose. */
/* Answer is returned in matrix B.                        */
void matmultLLT(real **L, real **B, int n)
{
  int i, j, k;

  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      B[i][j] = 0.0;
      for(k=0; k<=j; k++)
	B[i][j] += L[i][k] * L[j][k];

      if(j < i)
	B[j][i] = B[i][j];
    }
}

/* Use cholesky decomposition to get log determinant */
/* L should already be in Cholesky form */
real log_determinant(real **L, int n)
{
  real det;
  int i;

  det = log(L[0][0]);
  for(i=1; i<n; i++)
    det += log(L[i][i]);
  return(2.0*det);
}

/* return inner product of two vectors */
real dotprod(real *x, real *y, int n)
{
  int i;
  real z=0.0;

  for(i=0; i<n; i++)
    z += x[i]*y[i];

  return(z);
}
/****************************************************************
 * Solve
 *                          b = L * x
 *
 * where b and x are column vectors, and L is lower triangular with
 * diagonal entries not necessarily equal to 1.0.  The diagonal entries
 * of L *ARE* referenced.  On return b contains the solution vector x;
 ****************************************************************/
int col_L_solve(real **L, int n, real *b)
{
    int i;
    register int j;
    register real x;
    
    for(i=0; i<n; i++)
    {
      if(L[i][i] <= 0.0)
	return(1);
      b[i] = x = b[i]/L[i][i];
      for(j=i+1; j<n; j++)
	b[j] -= x * L[j][i];
    }
    return(0);
}
/*****************************************************************
 * Solve 
 *                            b' = x' * L
 * where b' and x' are row vectors, and L is lower triangular with diagonal
 * entries not necessarily equal to 1.0.  The diagonal entries of L *ARE*
 * referenced.  On return b contains the solution vector x;
 *****************************************************************/
int row_L_solve(real **L, int n, real *b)
{
    int i;
    register int j;
    register real x;
    
    for(i=n-1; i>=0; i--)
    {
      if(L[i][i] <= 0.0)
	return(1);
      b[i] = x = b[i]/L[i][i];
      for(j=0; j<i; j++)
	b[j] -= x * L[i][j];
    }
    return(0);
}

/****************************************************************
Mahalanobis distance: x'*inv(A)*x.
Should be used after a call to cholesky.
ON INPUT:
x	a vector of n reals
L	an nXn matrix of reals, as allocated by alloc2d.  The lower
        triangle of a should contain the left Cholesky factor of A.
n	dimension of square matrix L.
ON RETURN:
Function returns x'*inv(a)*x if successful.  This value is always positive,
	and negative returns indicate an error.
L	is unchanged.
x	is destroyed.
ERROR RETURN:
Returned value is < 0 if Cholesky factor is singular.

ALGORITHM:

The Cholesky decomposition turns
                  -1
        q = x' * A  * x.
 
into
                  -1   -1         2
        q = x' * L'  * L * x = |z|
where
             -1
        z = L * x
 
To find z we need only solve the triangular system
 
        L * z = x.
*****************************************************************/
real mahal(real *x, real **L, int n)
{
  int i;
  real dist=0.0;

  /* solve triangular system */
  if(col_L_solve(L, n, x))
  {
    fprintf(stderr,"\nmahal: cholesky factor is singular");
    exit(1);
  }

  /* get squared norm of solution vector */
  dist = 0.0;
  for(i = 0; i<n; i++)
    dist += x[i]*x[i] ;
  
  return(dist);
}

/***************************************************************

On input L is the left cholesky factor of a n X n matrix of
reals, called "a".  In other words, L is a lower triangular matrix
such that L * L' = a.   

On return, L contains the inverse of the input matrix L, and
s contains the inverse of "a", calculated as inv(L)' * inv(L).

Sources: Numerical Recipes, from p 312 of the Linpack Users' Guide,
1979 edition, and the Linpack routine spodi.f. 
***************************************************************/
void cholinv(real **L, real **s, int n)
{
  int i, j, k;
  real sum;

  /* invert L in place */
  for(i=0; i<n; i++)
  {
    L[i][i] = 1.0/L[i][i];
    for(j=i+1; j<n; j++)
    {
      sum = 0.0;
      for(k=i; k<j; k++)
	sum -= L[j][k] * L[k][i];
      L[j][i] = sum/L[j][j];
    }
  }

  /* form inv(A) = inv(L)' * inv(L) */
  for(i=0; i<n; i++)
  {
    for(j=0; j<=i; j++)
    {
      s[i][j] = 0.0;
      for(k=i; k<n; k++)
	s[i][j] += L[k][i] * L[k][j];
    }
    for(j=i+1; j<n; j++)
    {
      s[i][j] = 0.0;
      for(k=j; k<n; k++)
	s[i][j] += L[k][i] * L[k][j];
    }
  }
}
