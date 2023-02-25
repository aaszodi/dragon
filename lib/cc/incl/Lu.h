#ifndef __LU_HEADER__
#define __LU_HEADER__

// ==== HEADER Lu.h ====

/* LU-decomposition and linear equation solver
 * routines for square matrices.
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ---- STANDARD HEADERS ---- 

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

// ---- INCLUDE FILES ----

#include "Vector.h" 
#include "Sqmat.h"
#include "Array.h"

// ==== CLASSES ====

/* Lu_: a class for doing LU-decomposition. Can be asked to
 * calculate determinants and solve linear equation systems.
 */
class Lu_
{
    // data
    private:
    
    Sqmat_ Lu;	// holds the decomposition
    Array_<unsigned int> Perm;	// permutation vector
    int Psign;	// permutation sign
    
    // methods
    public:
    
	// constructors
    /* Inits to an N x N problem (default 2). */
    Lu_(unsigned int N=2): Lu(N), Perm(N), Psign(0) {}
    
	// decomposition
    /* decomp(): performs an LU-decomposition on the square matrix A.
     * Return value: the permutation sign (0 if A was singular).
     */
    int decomp(const Sqmat_& A);
    
    /* det: calculates the determinant from the decomposition results
     * stored inside. It is assumed that decomp() has been called beforehand.
     */
    double det() const;
    
	// solution
    /* solve(): solves the linear equation A*x=B and returns x.
     * The calling object is supposed to have been primed with
     * decomp(A) before this call. B is the "right-hand-side"
     * vector which is preserved. Checks for dim mismatches and returns
     * B untouched if not everything is OK.
     */
    Vector_ solve(const Vector_& B) const;
    
    /* lineq(): solves the nxn linear equation A*X=B and returns X.
     * Iterative improvement is done a la Num. Recipes (max. Maxit times, 
     * default 0).
     * This routine should be used if only one right-hand side vector B
     * is present. For several right-hand sides, use the decomposition
     * and solver routines separately. 
     * Return value: 1 means OK, 0 means A was singular, 
     * -1 means a dim mismatch occurred. 
     */
    int lineq(const Sqmat_& A, const Vector_& B, Vector_& X, unsigned int Maxit=0);
    
};
// END OF CLASS Lu_

// ==== END OF HEADER Lu.h ====

#endif
