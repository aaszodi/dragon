#ifndef MATBASE_CLASS
#define MATBASE_CLASS

// ==== HEADER Matbase.h ====

/* Abstract base class for rectangular and triangular matrices. */

/* MATRIX LIBRARY CLASS FAMILY TREE
 * 
 *         [ Matbase_ ]
 *               |
 *        +------+-------+
 *        |              |
 *        V              V
 *  [ Rectbase_ ]   [ Sqbase_ ]
 *        |              |
 *  +-----+-----+  +-----+-----+
 *  |           |  |           |
 *  V           V  V           V
 * Matrix_     Sqmat_      Trimat_
 */

// SGI C++ 4.0, IRIX. 5.3, 4. May 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>

// ---- INCLUDE FILES ----

#include "Vector.h"

// ==== CLASSES ==== 

/* Class Matbase_ : an abstract matrix class that serves as a base class 
 * for all other specialised matrix classes.
 * The elements of the matrix are stored in a one-dimensional double
 * array in row major order. A second array of double pointers, 
 * Rows, holds pointers to each row. The element array itself is pointed to
 * by Rows[0]. This arrangement can be used as the traditional 
 * array-of-pointers approach but here only two arrays have to be
 * allocated and copy operations can be performed on contiguous arrays.
 * Row numbers are stored in R (the column sizes are unknown at this
 * stage). Eno keeps track of the number of elements.
 */
class Matbase_
{
    // data
    protected:
    
    double **Rows;  // array of pointers to rows; Rows[0] points to the elements
    unsigned int R, Eno;	// row numbers, total no. of elements
    
    // public methods interface
    public:
    
	// constructors
    /* Init a matrix to have Row>=1 rows (default 3). The element array
     * is not allocated and Rows[] is not initialised.
     * The element number is set to 0.
     */
    Matbase_(unsigned int Row=3);
    
    /* Copy constructor: allocates the Rows array to the size of Mb,
     * does not initialise it, sets R to Mb.rno(), sets Eno=0.
     */
    Matbase_(const Matbase_& Mb);
    
	// destructor
    virtual ~Matbase_() { delete [] Rows; }
    
	// size
    unsigned int rno() const { return(R); }
    virtual unsigned int cno() const =0;
    
	// access (unsafe)
    const double *operator[](unsigned int Idx) const { return(Rows[Idx]); }
    double *operator[](unsigned int Idx) { return(Rows[Idx]); }
    
    /* Safe access: the overlaid function call operator is used for this.
     * Error messages are printed if one of the two indices is out of
     * range and the index in question will be replaced by 0.
     * Note that there are two versions of this function; one for
     * const access, the other for modifiable lvalues. Both are pure
     * virtual at this stage
     */
    virtual double operator()(unsigned int i, unsigned int j) const =0;
    virtual double& operator()(unsigned int i, unsigned int j) =0;
    
    /* row(Idx): returns the Idx-th row as a Vector_ object safely.
     * row(Vec,Idx): sets the elements in the Idx-th row to the elements
     * of the vector Vec, provided Idx is legal and Vec has the right
     * number of elements.
     */
    virtual Vector_ row(unsigned int Idx) const =0;
    virtual void row(const Vector_& Vec, unsigned int Idx) =0;
    
    /* col(Idx): returns the Idx-th column as a Vector_ object safely.
     * col(Vec,Idx): sets the elements in the Idx-th column to the elements
     * of the vector Vec, provided Idx is legal and Vec has the right
     * number of elements.
     */
    virtual Vector_ col(unsigned int Idx) const =0;
    virtual void col(const Vector_& Vec, unsigned int Idx) =0;
    
    /* set_values(): sets all elements to a specified value Val.
     * Default is 0.0.
     */
    void set_values(double Val=0.0);
    
    /* get_array(): converts the calling object to a traditional
     * double ** array and sets Row and Col to its corresponding
     * dimensions. Returns a NULL ptr if allocation failed.
     */
    virtual double** get_array(unsigned int& Row, unsigned int&Col) const =0;
    
	// FORTRAN indexing
    
    /* ftn_idx(): moves the data pointers so that 1..N
     * indexing could be used instead of 0..N-1. Note that
     * this operation is inherently unsafe (the "safe indexing"
     * method knows nothing about it!) and should be used only
     * when writing wrappers around Numerical Recipes code etc.
     * and only if the unsafe [][] indexing is used.
     * Make sure the object abused by this method is not destroyed
     * and indexing is reset via c_idx() as soon as possible.
     */
    void ftn_idx();
    
    /* c_idx(): resets the ptrs to their normal C-style position
     * after a call to ftn_idx().
     */
    void c_idx();
    
	// arithmetics
    /* Matrix addition and subtraction (+,+=,-,-=) check for dimension
     * mismatches. If a mismatch is found then the left
     * operand is returned unchanged. Division by scalar ( /,  /=)
     * is always checked for division by zero and the unchanged matrix
     * is returned if a div by 0 is attempted. These routines operate
     * on all elements of the element array, regardless of the row/col
     * layout and therefore are defined in the base class here and
     * will simply be inherited.
     */
    Matbase_& operator+=(const Matbase_& Mat);
    Matbase_& operator-=(const Matbase_& Mat);
    Matbase_& operator*=(double Factor);
    Matbase_& operator/=(double Div);
    
	// linear algebra
    /* Matrix*vector multiplication: the vector Vec is returned if a
     * dimension mismatch was detected.
     */
    virtual Vector_ operator*(const Vector_& Vec) const =0;

	// Printing
    /* list_matrix: lists calling object to stdout with entries occupying Width chars,
     * in scientific format, Prec digits precision. If a row takes up 
     * more than Linewidth chars, then the matrix is cut up nicely. 
     * If Width<Prec+8,  then it is adjusted but no warning is given.
     * Pure C++ I/O. See the redefined << in the Matbase.h header.
     * Could be called by all derived classes.
     */
    virtual void list_matrix(ostream& Out, 
	    unsigned int Prec=2, unsigned int Width=10, 
	    unsigned int Linewidth=80) const;
    
	// end of public interface
    
    // protected methods
    protected:
    
	// printing
    virtual void print_rows(ostream& Out, unsigned int Sizew, unsigned int Jbeg, 
	unsigned int Items, unsigned int Width, unsigned int Prec) const=0;
    
	// memory management
    double** alloc_rows(unsigned int Rowno);
    double* alloc_elems(unsigned int Elemno);
    virtual void init_rowptrs(unsigned int Colno) =0;
    
        // error messages
    enum Errtype_ { NO_MEM, DIV_BY_ZERO, BAD_ROWRANGE, BAD_COLRANGE, DIM_MISMATCH };
    void prt_err(const Errtype_ Etyp, const char *Funcnm) const;
    
};
// END OF CLASS Matbase_

// ---- GLOBAL PROTOTYPES ----

/* Overloaded version of <<: always calls the proper virtual
 * function for matrix listing.
 */
ostream& operator << (ostream& Ostr, const Matbase_& Mat);

// ==== END OF HEADER Matbase.h ====

#endif	/* MATBASE_CLASS */
