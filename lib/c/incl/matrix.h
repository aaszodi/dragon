#ifndef MATRIX_HEADER
#define MATRIX_HEADER

/* ==== HEADER matrix.h ==== */

/* Header for square and lower triangle matrices: the latter
 are stored economically. */

/* ANSI C, IRIX 5.2, 5. Aug. 1994. Andris Aszodi */

/* ---- HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <limits.h>	/* for precision stuff */
#include <float.h>
#include <math.h>

/* ---- GLOBAL TYPES ---- */

/* WARNING: The routines do not check whether the matrix is triangular
 or square. The typedef's below might help during compilation. */

typedef double **Matrix_ ;	/* general type */
typedef double **Trimat_ ;	/* triangular */
typedef double **Sqmat_ ;	/* square */

/* ---- PROTOTYPES ---- */

/* alloc_trimat: allocates space for a triangular matrix with Size rows.
 * The triangle contains the main diagonal as well. 
 * Returns the pointer to the matrix or NULL if alloc failed.
 */
Trimat_ alloc_trimat(unsigned int Size);

/* free_matrix: frees up space occupied by Mat. Use this routine for
 * both triangular and square matrices. Mat itself is not changed.
 */
void free_matrix(double **Mat);

/* list_trimat: lists Mat to stdout with entries occupying Width chars,
 * Prec digits precision. If a row takes up more than Linewidth chars,
 * then the matrix is cut up nicely.
 */
void list_trimat(Trimat_ Mat, int Size, int Linewidth,
	int Width, int Prec);

/* alloc_sqmat: allocates space for a square matrix (Size*Size).
 * Returns the pointer to the matrix or NULL if alloc failed.
 */
Sqmat_ alloc_sqmat(unsigned int Size);

/* list_sqmat: lists Mat to stdout with entries occupying Width chars,
 * Prec digits precision. If a row takes up more than Linewidth chars,
 * then the matrix is cut up nicely.
 */
void list_sqmat(Sqmat_ Mat, int Size, int Linewidth,
	int Width, int Prec);

/* lu_decomp: performs an LU-decomposition in place on the n*n matrix
 * A. Based on partial pivoting: the row permutations are done in Perm[]
 * (assumed to be of correct size) and will be used by lu_solve().
 * If Perm==NULL, then a permutation vector will be used
 * internally and will be freed before return. This option can be
 * used when only the determinant is calculated from the LU-decomposition.
 * Return value: the sign of the determinant of the permutation
 * matrix (+/-1) or 0 if A is singular or n<=0.
 */
int lu_decomp(Sqmat_ A, int n, int *Perm);

/* lu_det: calculates the determinant of the n x n LU-decomposed
 * square matrix Lu. Psign is the permutation sign returned by
 * lu_decomp().
 */
double lu_det(const Sqmat_ Lu, int Psign, int n);

/* lu_solve: solves the linear equation A*x=b (A is n*n, x,b are n long).
 * A is supposed to have been LU-decomposed by lu_decomp() above and
 * the row permutation is stored in Perm[]. b[] is the "right-hand-side"
 * vector which contains the solution on return.
 */
void lu_solve(const Sqmat_ A, const int Perm[], double b[], int n);

/* ==== END OF HEADER matrix.h ==== */

#endif		/* MATRIX_HEADER */
