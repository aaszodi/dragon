#ifndef SQMAT_CLASS
#define SQMAT_CLASS

// ==== HEADER Sqmat.h ====

/* Double-precision square matrix class. Derived from the
 * general RxC rectangular matrix class Rectbase_ and the
 * square matrix base class Sqbase_ . 
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris 

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

// ---- INCLUDE FILES ----

#include "Sqbase.h"	// the abstract base classes
#include "Rectbase.h"
#include "Matrix.h"	// the RxC matrix class (for conversions)
#include "Vector.h"	// vectors

// ==== CLASSES ====

/* Class Sqmat_ : double-precision square matrices. Derived
 * from the general RxC rectangular matrix class. Storage
 * is exactly the same as well as most methods. A few operations
 * are added such as diagonal manipulations, trace and in-place
 * transposition. The main reason for defining a separate 
 * square matrix class is to implement linear equations, 
 * determinants, eigenroutines etc. 
 */
class Sqmat_ : public Rectbase_, public Sqbase_
{
    // methods interface
    public:
    
	// constructors
    /* Allocates a Size x Size matrix (default size is 3x3).
     * Initialises all elements to 0.0.
     */
    Sqmat_(unsigned int Size=3) : Matbase_(Size), Rectbase_(Size) {}
    
    /* Allocates a square matrix to hold the elements of the
     * traditional array Arr. It is the caller's responsibility
     * to make sure that Arr is indeed Size x Size. If Arr==NULL
     * then the matrix elements will be set to 0.0
     */
    Sqmat_(const double **Arr, unsigned int Size) : Matbase_(Size), Rectbase_(Arr, Size, Size) {}
    
    /* Copy constructor */
    Sqmat_(const Sqmat_& Sq): Matbase_(Sq), Rectbase_(Sq) {}
    
    /* Init with a Rectbase_ (RxC) object (conversion ctor). The resulting Sqmat_
     * will contain the full matrix so that its size will be max(R, C)
     * with the extra rows/cols padded with zeroes. Also works with Sqmat_.
     */
    Sqmat_(const Rectbase_& Rbase);
    
	// assignment
    Sqmat_& operator=(const Sqmat_& Sq);
    
	// size
    /* set_size(): resets the size to Size x Size. The upper left
     * corner overlap is preserved, if the new size is
     * less than the original then the extra rows and cols will be lost.
     * If new rows/cols are added then they will be initialised to 0.0.
     * Zero row/col numbers not permitted.
     * If the new size is 0 or the same as the old then no action is taken.
     */
    void set_size(unsigned int Size);
    
	// linear algebra
    /* Matrix multiplication: the left operand is returned if
     * a dimension mismatch is detected. The code is exactly the same
     * as the corresponding Matrix_ operator; however, the return type
     * does not allow us to decl/def the whole thing in Rectbase_
     * (where it logically belongs).
     */
    Matrix_ operator*(const Rectbase_& Mat) const;
    
    /* Matrix*vector multiplication: the vector is returned if a
     * dimension mismatch was detected.
     */
    Vector_ operator*(const Vector_& Vec) const;
    
	// arithmetics
    /* Matrix addition and subtraction (+,-,+=,-=) check for dimension
     * mismatches. If a mismatch is found then the left
     * operand is returned unchanged. Division by scalar ( / )
     * is always checked for division by zero and the unchanged matrix
     * is returned if a div by 0 is attempted. *=, /= are inherited
     * from Matbase_ (could not be declared there because of the
     * return type).
     */
    Sqmat_ operator+(const Sqmat_& Mat) const;
    Sqmat_ operator-(const Sqmat_& Mat) const;
    Sqmat_ operator*(double Factor) const;
    friend Sqmat_ operator*(double Factor, const Sqmat_& Mat);
    Sqmat_ operator/(double Div) const;
    
	// transpose
    /* get_transpose(): returns the transposed calling object. */
    Sqmat_ get_transpose() const;
    
    /* transpose_inplace(): transposes the calling object "in place". */
    void transpose_inplace();

};
// END OF CLASS Sqmat_

// ==== END OF HEADER Sqmat.h ====

#endif	/* SQMAT_CLASS */
