#ifndef _NUMERICS_H_
#define _NUMERICS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>


#include "library.h"

#define NRMAX(i,j) ( (i)<(j) ? (j):(i) )
#define NRMIN(i,j) ( (i)<(j) ? (i):(j) )
#define NRTINY 1.0e-10


// **********************************************
//  float **  basic functions
// **********************************************

float ** allocate_float_matrix(int nrows, int ncols);

void desallocate_float_matrix(float **matrix, int nrows, int ncols);

void print_float_matrix(float **matrix, int nrows, int ncols);

void product_square_float_matrixes(float **result,float **matrix1,float **matrix2,int n);

void float_vector_matrix_product(float **a, float *x,float *y,int n);

void order_decreasing(float *values, int *indexos, int size);

// **********************************************
//  LU based algorithms
// **********************************************

// Solves Ax=b by using lu decomposition 
// a matrix a[1..n][1..n] is replaced by the LU decompositions of a rowwise permutation of itself
// b[1..n] and x[1..n]
int lusolve(float **a, float *x, float *b, int n);


// Computes the inverse of a column by column using the LU decomposition
// a matrix a[1..n][1..n] is replaced by the LU decompositions of a rowwise permutation of itself
int luinv(float **a, float **inv, int n);


/*-- LU decomposition,  From Numerical Recipes in C */
/* Given a matrix a[1..n][1..n] this routine replacess it by the LU decompositions of a rowwise permutation of itself. */
/* a and n are input, a is output, arranged as in equation (2.3.14) above; indx[1..n] in an output vector that records */
/* the row permutation effected by the partial pivoting; d is output as +-1 depending on whether the number of row     */
/* interchanges was even or odd respectively.                                                                          */
int ludcmp(float **a, int n, int *indx, float *d); /* LU decomposition */

/* Solves the set of n linear equations Ax=b. Here a[0..n-1][0..n-1] as input, not as the matrix A but rather as its LU decomposition,*/
/* determined by the routine ludcmp. indx[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is input as the    */
/* right hand side vector and returns with the solution vector x. */
void lubksb(float **a, int n, int *indx, float *b); /* LU linear solution */



// **********************************************
//  Householder reduction to tridiagonal matrix
// **********************************************

/*- Householder reduction of a real symmetric matrix A[0..n-1][0..n_1]. On output, A is ---*/
/*- replaced by the ortogonal matrix Q effecting the transformation. d[0..n-1] returns ----*/
/*- the diagonal elements of the diagonal matrix, and e[0..n-1] the off-diagonal elements -*/
/*- with e[0]=0. 									   */
void symmetric_2_tridiag_HH(float **a,float *d,float *e,int n);



// **********************************************
//  QL decomposition algorithm to be used with tridiagonal matrices
// **********************************************

/* QL with implicit shifts, to determine eigenvalues and eigenvectors of  a real, symmetric, tridiagonal matrix */
/* d[0..n-1] contains the diagonal elements of the matrix and as output the returns the eigenvalues             */
/* e[0..n-1] contains the sub-diagonal elements with e[0] arbitrary, on output e is destroyed			*/
/* z[0..n-1][0..n-1] contain the identity matrix in input or the output of symmetric_2_tridiag_HH if previously */
/* applied. 													*/

int eigenvalues_QLI(float * d, float *e,float **z, int n);


// **********************************************
//  Singular Value Decomposition
// **********************************************

//  m_U (A.nrow(), A.ncol)(), m_V(A.ncol(),A.ncol()), m_W(A.ncol())
void compute_svd(float **A, float **m_U, float **m_V, float *m_W, int rows, int cols);


// **********************************************
//  PCA
// **********************************************


void  pca_center_data(float **X,float *baricenter,int n, int p);


////////////////////// SVD
// Data is already centered by using pca_center_data
// X:  n rows and p columns
// p:  dimension
// n:  number of samples
// S:  float[p] returns the eigenvalues
// V:  float[p][p] returns the principal vectors
// U:  float[n][p] returns the new coefficients
void compute_pca_svd(float **X,float *S,float **V, float **U,int n, int p);



///////////////////// FORCE BRUTE
int  compute_pca(float **X,float **Pcs, float *sVar, int n, int p);
void compute_coefic(float **X,float **PCs,int nvec,float **Coef, int n, int p);
float **  covariance_matrix(float **x,int n, int p);

// Data is already centered by using pca_center_data
// X:  n rows and p columns
// p:  dimension
// n:  number of samples
// S:  float[p] returns the eigenvalues
// V:  float[p][p] returns the principal vectors
// U:  float[n][p] returns the new coefficients
void compute_pca_brute(float **X,float *S,float **V, float **U,int n, int p);

#endif

