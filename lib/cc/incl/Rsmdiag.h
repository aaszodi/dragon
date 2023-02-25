#ifndef RSMDIAG_CLASS
#define RSMDIAG_CLASS

// ==== PROJECT DRAGON: HEADER Rsmdiag.h ====

/* Real symmetric matrix diagonalisation class.
 * Use if only selected eigenvectors are wanted.
 * Based on NETLIB routines.
 */

// SGI C++ 7.1, IRIX 6.2, 21. June 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <math.h>

// ---- MODULES ----

#include "Sqmat.h"
#include "Trimat.h"
#include "Vector.h"
#include "Safety.h"

/* NOTE: The SGI N32/N64 compilers recognise the built-in bool type
 * and define the macro _BOOL. This workaround is provided for
 * the O32 compiler.
 */
#if !(defined(_BOOL) || defined(HAS_SEEN_BOOL))
#define HAS_SEEN_BOOL
typedef unsigned char bool;
static const bool false=0;
static const bool true=1;
#endif	/* _BOOL */

// ==== CLASSES ====

/* Rsmdiag_: a class implementing symmetric matrix diagonalisation.
 * All eigenvalues are generated in decreasing order, and then the
 * user can request the first m eigenvectors. Note that this is
 * ineffective unless m<Size/4.
 */
class Rsmdiag_
{
    // constants
    private:
    static const Safety_ Safe; // division and hypot() safeguard
    
    // data
    private:
    
    Trimat_ Qmat;   // work matrix
    double *d, *e, *e2, *w;	// tridiagonal matrix storage and eigenvalue temp
    int *Index;	    // submatrix index array
    bool Ftnidx;   // true if indices are shifted
    
    // methods
    public:
    
	// ctor, dtor
    Rsmdiag_(): 
	Qmat(), d(NULL), e(NULL), e2(NULL), w(NULL), 
	Index(NULL), Ftnidx(false) {}
    
    ~Rsmdiag_() { set_size(0); }
    
	// diagonalisation
    /* get_evals(): obtain all eigenvalues of a symmetric matrix Mat
     * and put them into the Evals vector (size set if necessary).
     * The eigenvalues are sorted in decreasing order.
     * Return value: 0 if OK, k>0 if the k:th eigenvalue failed to
     * converge.
     */
    int get_evals(const Trimat_& Mat, Vector_& Evals);
    
    /* get_evecs(): obtain the first Evno eigenvectors. It is assumed but
     * not checked that a corresponding get_evals() has been executed
     * beforehand. The eigenvectors will be placed into the first Evno
     * columns of the Evecs matrix whose size will be set silently.
     * Return value: 0 if OK, r>0 if the r:th eigenvector failed
     * to converge.
     */
    int get_evecs(int Evno, Sqmat_& Evecs);
    
    // hidden methods
    protected:
    
	// diag from NETLIB
    void tred_1();
    int imt_qlv();
    int inv_iter(int m, Sqmat_& z) const;
    void tr_bak1(int m, Sqmat_& z) const;

	// auxiliaries
    static double epsilon(double x);
    void c_idx();
    void ftn_idx();
    unsigned int set_size(unsigned int Size);
    
    // forbidden methods
    private:
    Rsmdiag_(const Rsmdiag_&);	// no copy ctor
    Rsmdiag_& operator=(const Rsmdiag_&);   // no assignment
    
};
// END OF CLASS Rsmdiag_

// ==== END OF HEADER Rsmdiag.h ====

#endif	/* RSMDIAG_CLASS */
