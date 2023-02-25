#ifndef __MATRIX_CLASS__
#define __MATRIX_CLASS__

// ==== HEADER Matrix.h ====

/* Double-precision class for general RxC rectangular matrices. 
 * The Vector_ class is used for various linear algebra operations.
 */
 
// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

// ---- INCLUDE FILES ---- 

#include "Matbase.h"
#include "Rectbase.h"

// ==== CLASSES ====

/* Class Matrix_ : a general RxC rectangular matrix class.
 * The elements of the matrix are stored in a one-dimensional double
 * array in row major order. A second array of double pointers, 
 * Rows, holds pointers to each row. The element array itself is pointed to
 * by Rows[0]. Derived from Rectbase_.
 */
class Matrix_ : public Rectbase_
{
    // data
    protected:
    unsigned int C;	// number of columns
    
    // methods interface
    public:
    
	// constructors
    /* Allocates a Row x Col matrix. Without any arguments, a 3x3
     * matrix will be created. With one argument, the Col parameter
     * defaults to 0,  which is interpreted as an instruction to create
     * a square Row x Row matrix. Square matrices may be created
     * by Matrix_(N, N) or Matrix_(N) or Matrix_(N, 0).
     * All elements will be initialised to 0.0.
     */
    Matrix_(unsigned int Row=3, unsigned int Col=0): Matbase_(Row), Rectbase_(Row, Col) { C=Eno/R; }
    
    /* Allocate a Row*Col matrix as above and initialise it with
     * a traditional double ** style array Arr. It is the caller's 
     * responsibility to make sure that Arr has the correct dimensions.
     * If Arr==NULL then the matrix elements will be set to 0.0
     */
    Matrix_(const double **Arr, unsigned int Row, unsigned int Col):
	Matbase_(Row), Rectbase_(Arr, Row, Col), C(Col) {}
    
    /* Init with a base ref. (base->derived conversion ctor) */
    Matrix_(const Rectbase_& Rbase): Matbase_(Rbase.rno()), Rectbase_(Rbase) { C=Eno/R; }
    
    /* the copy constructor */
    Matrix_(const Matrix_& Mat): Matbase_(Mat), Rectbase_(Mat), C(Mat.cno()) {}
    
    /* Assignment */
    Matrix_& operator=(const Matrix_& Mat);
    
	// size
    /* set_size(): resets the size to Rno x Cno. The upper left
     * corner overlap is preserved, if the new row/col number is
     * less than the original then the extra rows/cols will be lost.
     * If new rows/cols are added then they will be initialised to 0.0.
     * Zero row/col numbers not permitted.
     * If the new size is 0 or the same as the old then no action is taken.
     */
    void set_size(unsigned int Rno, unsigned int Cno);
    unsigned int cno() const { return(C); }
    
	// simple linear algebra
    /* Matrix multiplication: the left operand is returned if
     * a dimension mismatch is detected. There is no in-place (*=) version.
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
    Matrix_ operator+(const Rectbase_& Mat) const;
    Matrix_ operator-(const Rectbase_& Mat) const;
    Matrix_ operator*(double Factor) const;
    friend Matrix_ operator*(double Factor, const Rectbase_& Mat);
    Matrix_ operator/(double Div) const;
    
	// transposition
    /* get_transpose(): returns the transpose of the calling object. */
    Matrix_ get_transpose() const;
    
};
// END OF CLASS Matrix_

// ==== END OF HEADER Matrix.h ====

#endif
