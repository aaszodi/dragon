#ifndef SVD_CLASS
#define SVD_CLASS

// ==== HEADER Svd.h ====

/* Singular value decomposition based on the algorithm
 * in Numerical Recipes.
 * The two matrices and the weight vector produced by SVD
 * are bundled into a little class.
 */

// SGI C++, IRIX 6.2, 20. June 1998. Andris

// ----	STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <math.h>

// ---- INCLUDE FILES ----

#include "Matrix.h"
#include "Sqmat.h"
#include "Trimat.h"
#include "Vector.h"
#include "Vmutils.h"
#include "Safety.h"

// ---- DEFINITIONS ----

#ifdef FLT_EPSILON
#define SVD_EPSILON (10.0*FLT_EPSILON)
#else
#define SVD_EPSILON (1e-7)
#endif

// ==== CLASSES ====

/* Class Svd_ : a simple class for Singular Value Decomposition.
 * SVD decomposes a rectangular matrix A into the form UWV'
 * where W is diagonal and V is square. U and V are Matrix_
 * and Sqmat_ objects, respectively, W is stored as a Vector_
 * object. The class provides a few methods for performing SVD and analysing
 * the results. A general linear equation solver is also provided.
 */
class Svd_
{
    // data
    private:
    
    static const Safety_ Safe;	// safe division and hypot() classlet
    
    Matrix_ U;	// general Row x Col matrix
    Vector_ W;	// Col-long vector
    Sqmat_ V;	// Col x Col square matrix
    unsigned int R, Rorig, C;	// row and column nos.
    
    // methods
    public:
    
	// constructors
    /* Set up SVD for a Row x Col matrix where Row>=Col. If Row<Col,
     * then the rows of U will be padded to make it ColxCol
     * with a warning. If either Row or Col is 0, will be changed
     * to 3 (cf matrix/vector ctors) and a warning printed. 
     * The default size is also 3x3.
     * The members will be initialised to 0.0 (all elements).
     */
    Svd_(unsigned int Row=3, unsigned int Col=3);
    
	// access
    /* u(), w() and v() return const refs to the corresponding components. */
    const Matrix_& u() const { return(U); }
    const Vector_& w() const { return(W); }
    const Sqmat_& v() const { return(V); }
    
    /* rno() and cno() return the current row and col. numbers */
    unsigned int rno() const { return(R); }
    unsigned int cno() const { return(C); }
    
    /* set_size(): modifies the sizes of the SVD components to
     * accommodate a Row x Col matrix. If Row < Col, then the
     * rows of U will be padded to give a Col x Col matrix.
     * If Row or Col is 0, then 3 will be used instead and a warning printed.
     * Without arguments, the size is set to 3x3.
     */
    void set_size(unsigned int Row=3, unsigned int Col=3);
    
	// SVD routines
    /* make_decomp(): the SVD according to Num. Recipes, done
     * on the general matrix A (which is preserved). The
     * three member objects U, W, V will be set in the calling
     * object. Return value: 1 if the iteration limit (hard-coded
     * in eigen_ql()) is exceeded, 0 if OK, -1 if a dimension
     * mismatch occurred.
     */
    int make_decomp(const Matrix_& A);
	
    /* rank_cond(): checks the N singular values W[] of a matrix 
     * after SVD. The condition number in *Cond
     * (ratio of the largest and smallest singular value) is also
     * calculated if Cond!=NULL. The singular values which are smaller than
     * Eps times the largest are set to 0.0.
     * Return value: the rank of the matrix.
     */
    unsigned int rank_cond(double Eps=SVD_EPSILON, double *Cond=NULL);
    
    /* lin_solve(): back-substitution routine for solving linear equations
     * AX=B. A should be SV-decomposed into U, W and V' by make_decomp()
     * and the weight vector should be "conditioned" (small entries
     * zeroed) by rank_cond() prior to the call to this routine.
     * Return value: the X vector.
     */
    Vector_ lin_solve(const Vector_& B) const;
    
    /* reset_data(): zeroes the matrices and the W vector. */
    void reset_data() { U.set_values(); W.set_values(); V.set_values(); }
    
    /* <<: the overloaded output operator. Prints a neat listing to Out */
    friend ostream& operator<<(ostream& Out, const Svd_& Svd);
    
    private:
    int svd_core();
    void utb(const Vector_& B, Vector_& Utb) const;

    static double safe_div(double Num, double Denom, int Lineno=0);
    
};
// END OF CLASS Svd_

// ==== END OF HEADER Svd.h ====
#endif	/* SVD_CLASS */
