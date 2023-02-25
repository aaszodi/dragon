#ifndef QL_HEADER
#define QL_HEADER

// ====	HEADER Ql.h ====

/* Diagonalisation of symmetric real matrices (Trimat_ class).
 * Based on Numerical Recipes, with some modifications and
 * C++ -isation. Implements the QL-algorithm with implicit shifts.
 * Also contains a simple iteration by "popular demand".
 */

// SGI C++ 4.0, IRIX 5.3, 2. May 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream.h>

// ---- INCLUDE FILES ----

#include "Trimat.h"
#include "Sqmat.h"
#include "Vector.h"

/* ---- PROTOTYPES ---- */

/* eigen_ql(): a 'wrap' function driving the Housholder and QL routines.
 * Takes a lower triangular matrix Mat as input (preserved) and
 * produces the eigenvalues in Eval and the eigenvectors in Evec.
 * The sizes of these are adjusted silently if necessary.
 * Eigenvalues are sorted in decreasing order and the corresponding
 * eigenvectors are the _col_vectors_ of Evec.
 * Index shifts are performed to hack around the Fortranese
 * [1..N] convention of Numerical Recipes.
 * Return value: 0 if OK, 1 if iteration limit was exceeded in tqli().
 */ 
int eigen_ql(const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec);

/* eigen_positer(): tries to find the Poseno largest positive eigenvalues
 * and corresponding eigenvectors of Mat. Mat, Eval, Evec are
 * as in eigen_ql(). Poseno<=size of Mat. The actual number
 * of positive eigenvalues are returned.
 */
int eigen_positer(int Poseno, const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec);

/* eigen_poscheb(): generates the first Poseno eigenvalues and
 * eigenvectors. Syntax is the same as that of eigen_positer()
 * but the algorithm used is the Chebyshev iteration (more robust).
 */
int eigen_poscheb(int Poseno, const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec);

// ==== END OF HEADER Ql.h ====

#endif	/* QL_HEADER */
