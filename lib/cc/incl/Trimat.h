#ifndef __TRIMAT_CLASS__
#define __TRIMAT_CLASS__

// ==== HEADER Trimat.h ====

/* Triangular matrix class. Symmetric square matrices are
 * stored as lower triangles (main diagonal included)
 * in this class to save memory. Derived from Sqbase_. 
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

// ---- INCLUDE FILES ----

#include "Sqbase.h"  // base class 
#include "Sqmat.h"  // for conversions
#include "Matrix.h"

// ==== CLASSES ====

/* Class Trimat_ : this class implements a less memory-hungry
 * symmetric square matrix class. Only the main diagonal and the
 * lower triangle are stored. The element and row pointer storage
 * has the same layout as in the ultimate ancestor (Matbase), 
 * but the row pointers point to increasingly longer rows.
 */
class Trimat_ : public Sqbase_
{
    /* public methods interface */
    public:
    
	// Constructors
    /* Init to SizexSize (default 3x3), fill up with 0. */
    Trimat_(unsigned int Size=3);
    
    /* Allocate a Row*Row triangular matrix and initialise it with
     * a traditional double ** style array Arr. It is the caller's 
     * responsibility to make sure that Arr has the correct dimensions.
     * Arr is supposed to have the "triangular" layout already but
     * there's no way of checking it. 
     * If Arr==NULL then the matrix elements will be set to 0.0
     */
    Trimat_(const double **Arr, unsigned int Row);
    
    /* The copy constructor */
    Trimat_(const Trimat_& Tri);
    
    /* Tri<--Sq conversion constructor. Inits with Sq's lower triangle */
    Trimat_(const Sqmat_& Sq);
    
	// destructor
    virtual ~Trimat_() { delete [] Rows[0]; }
    
	// assignment
    Trimat_& operator=(const Trimat_& Tri);
    
	// Conversion
    /* Converts a Trimat_ object into a Sqmat_ object, i.e.
     * makes a full symmetric matrix out of the sparse representation.
     */
    operator Sqmat_() const;
	
    /* get_array(): converts the calling object to a traditional
     * double ** array and sets Row and Col to its corresponding
     * dimensions. Returns a NULL ptr if allocation failed.
     * The returned array will retain the "triangular" layout
     * (similar to my earlier C "Trimat_" type).
     */
    double** get_array(unsigned int& Row, unsigned int& Col) const;
    
	// Size
    /* set_size(): Sets the size to Size x Size. If the new size is
     * larger than the old then the new rows will be set to 0.0, 
     * if less then the upper triangle will be preserved and the extra
     * rows will be lost. Zero size not allowed (no action).
     */
    void set_size(unsigned int Size);
    
    /* Safe access: the overlaid function call operator is used for this.
     * Error messages are printed if one of the two indices is out of
     * range and the index in question will be replaced by 0.
     * If the column index is higher than the row index (j>i)
     * then the two will be swapped (Trimat_s are symm square matrices).
     * Note that there are two versions of this function; one for
     * const access, the other for modifiable lvalues.
     */
    double operator()(unsigned int i, unsigned int j) const;
    double& operator()(unsigned int i, unsigned int j);
    
    /* row(Idx): returns the Idx-th row as a Vector_ object safely.
     * row(Vec,Idx): sets the elements in the Idx-th row to the elements
     * of the vector Vec, provided Idx is legal and Vec has the right
     * number of elements.
     */
    Vector_ row(unsigned int Idx) const;
    void row(const Vector_& Vec, unsigned int Idx);
    
    /* col(Idx): returns the Idx-th column as a Vector_ object safely.
     * col(Vec,Idx): sets the elements in the Idx-th column to the elements
     * of the vector Vec, provided Idx is legal and Vec has the right
     * number of elements.
     */
    Vector_ col(unsigned int Idx) const;
    void col(const Vector_& Vec, unsigned int Idx);
    
	// simple linear algebra
    /* Matrix multiplication: the left operand is returned if
     * a dimension mismatch is detected. For simplicity, the
     * "safe indexing" is used here with the associated performance
     * penalty.
     * Returns a general RxC Matrix_.
     */
    Matrix_ operator*(const Matbase_& Mat) const;
    
    /* Matrix*vector multiplication: the vector Vec is returned if a
     * dimension mismatch was detected.
     */
    Vector_ operator*(const Vector_& Vec) const;

	// Arithmetics
    /* Matrix addition and subtraction (+,-) check for dimension
     * mismatches. If a mismatch is found then the left
     * operand is returned unchanged. Division by scalar ( / )
     * is always checked for division by zero and the unchanged matrix
     * is returned if a div by 0 is attempted. These routines are
     * based on the "in-place" operators (+=, -=, *=, /=) inherited
     * from Matbase_ (could not be declared there because of the
     * return type).
     */
    Trimat_ operator+(const Trimat_& Mat) const;
    Trimat_ operator-(const Trimat_& Mat) const;
    Trimat_ operator*(double Factor) const;
    friend Trimat_ operator*(double Factor, const Trimat_& Tri);
    Trimat_ operator/(double Div) const;

	// transpositions
    /* Since the transpose of a symmetric matrix is always equal
     * to itself, the following methods are very simple.
     */
    void transpose_inplace() { return; }
    Trimat_ get_transpose() const { return(*this); }

    /* protected methods */
    protected:
	// printing
    void print_rows(ostream& Out, unsigned int Sizew, unsigned int Jbeg, 
	unsigned int Items, unsigned int Width, unsigned int Prec) const;
	// Memory management
    void init_rowptrs(unsigned int Colno);
};
// END OF CLASS Trimat_

// ==== END OF HEADER Trimat.h ====
#endif
