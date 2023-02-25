#ifndef __SIVA_CLASS__
#define __SIVA_CLASS__

// ==== HEADER Siva.h ====

/* Singular value decomposition based on the algorithm
 * in Pal Rozsa's book (Linearis algebra es alkalmazasai).
 * The two matrices and the weight vector produced by SVD
 * are bundled into a little class but access to them is
 * public.
 */

// SGI C++ 3.2.1, IRIX 5.2, 21. Nov. 1994. Andris

// ----	STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <math.h>

// ---- INCLUDE FILES ----

#include "Matrix.h"
#include "Sqmat.h"
#include "Trimat.h"
#include "Vector.h"
#include "Ql.h"
#include "Vmutils.h"

// ---- DEFINITIONS ----

#ifdef FLT_EPSILON
#define SIVA_EPSILON (10.0*FLT_EPSILON)
#else
#define SIVA_EPSILON (1e-7)
#endif

// ==== CLASSES ====

/* Class Siva_ : a simple class for SIngular VAlue decomposition.
 * SVD decomposes a rectangular matrix A into the form UWV'
 * where W is diagonal and V is square. U and V are Matrix_
 * and Sqmat_ objects, respectively, W is stored as a Vector_
 * object. For the sake of simplicity, these data member objects
 * are "public" which goes against encapsulation. The class
 * provides a few methods for performing SVD and analysing
 * the results. A general linear equation solver is also provided.
 */
class Siva_
{
    // data
    public:
    
    Matrix_ U;	// general Row x Col matrix
    Vector_ W;	// Col-long vector
    Sqmat_ V;	// Col x Col square matrix
    unsigned int R, Rorig, C;	// row and column nos.
    
    // methods
    public:
    
	// constructors
    /* Set up SVD for a Row x Col matrix where Row>=Col. If Row<Col,
     * then the rows of U will be padded to make it ColxCol and a
     * warning is printed. If either Row or Col is 0, will be changed
     * to 1 (cf matrix/vector ctors) and another warning printed. 
     * The members will be initialised to 0.0 (all elements).
     */
    Siva_(unsigned int Row, unsigned int Col);
    
	// SVD routines
	
    /* make_decomp(): the SVD according to P. Rozsa, done
     * on the general matrix A (which is preserved). The
     * three member objects U, W, V will be set in the calling
     * object. Return value: 1 if the iteration limit (hard-coded
     * in eigen_ql()) is exceeded, 0 if OK, -1 if a dimension
     * mismatch occurred.
     */
    int make_decomp(const Matrix_& A);
	
    /* rank_cond(): checks the N singular values W[] of a matrix 
     * after SVD. The condition number Cond
     * (ratio of the largest and smallest singular value) is also
     * calculated. The singular values which are smaller than
     * Eps times the largest are set to 0.0.
     * Return value: the rank of the matrix.
     */
    unsigned int rank_cond(double Eps, double& Cond);
    
    /* lin_solve(): back-substitution routine for solving linear equations
     * AX=B. A should be SV-decomposed into U, W and V' by make_decomp()
     * and the weight vector should be "conditioned" (small entries
     * zeroed) by rank_cond() prior to the call to this routine.
     * Return value: the X vector.
     */
    Vector_ lin_solve(const Vector_& B) const;
    
    /* reset_data(): zeroes the matrices and the W vector. */
    void reset_data() { U.set_values(); W.set_values(); V.set_values(); }
};
// END OF CLASS Siva_

// ---- GLOBAL PROTOTYPES ----

/* <<: the overloaded output operator. Prints a neat listing to Out */
ostream& operator<<(ostream& Out, const Siva_& Svd);

// ==== END OF HEADER Siva.h ====

#endif
