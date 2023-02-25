#ifndef __VMUTILS_HEADER__
#define __VMUTILS_HEADER__

// ==== HEADER Vmutils.h ====

/* Various utilities for vectors and matrices which are
 * not member functions of the vector and matrix classes
 * but nevertheless deserve to be included in the VM library.
 */

// SGI C++ 3.2.1, IRIX 5.2, 21. Nov. 1994. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ---- INCLUDE FILES ----

#include "Vector.h"
#include "Matbase.h"
#include "Matrix.h"
#include "Sqmat.h"
#include "Trimat.h"

// ---- PROTOTYPES ----

    // transpose products
    
/* trans_prod(): calculates the product X'X of a matrix X.
 * The result is a Trimat_ because it is symmetric. Safe
 * indexing is used which may incur a slight performance
 * penalty but it is necessary for generality here as X
 * is a base class ref which may be either a rectangular
 * Matrix_ or a Trimat_. See below for specialised versions.
 */
Trimat_ trans_prod(const Matbase_& X);

/* trans_mprod(): calculates the product X'X of a rectangular Matrix_ X.
 * Result is a symmetric Trimat_. This routine should be slightly
 * faster for rectangulars than trans_prod() above but it cannot
 * be used for Trimat_s.
 */
Trimat_ trans_mprod(const Matrix_& X);

/* trans_wprod(): Calculates the product X'WX where X is a general
 * matrix and W is a diagonal matrix (represented by a Vector_).
 * Safe indexing is used (cf comments on trans_prod()). In case
 * of a dim mismatch between X and W then X'X is calculated.
 */
Trimat_ trans_wprod(const Matbase_& X, const Vector_& W);

/* trans_mwprod(): calculates the product X'WX of a rectangular matrix
 * X and a diagonal matrix represented by the Vector_ W. Less general
 * but slightly faster than trans_wprod() above.
 */
Trimat_ trans_mwprod(const Matrix_& X, const Vector_& W);

// ==== END OF HEADER Vmutils.h ==== 

#endif

