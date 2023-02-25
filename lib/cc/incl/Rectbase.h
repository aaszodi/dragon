#ifndef __RECTBASE_CLASS__
#define __RECTBASE_CLASS__

// ==== HEADER Rectbase.h ====

/* Abstract base class for rectangular matrices (square and non-square). */
 
// SGI C++ 4.0, IRIX 5.3, 4. May 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

// ---- INCLUDE FILES ---- 

#include "Matbase.h"
#include "Vector.h"

// ==== CLASSES ====

/* Class Rectbase_ : ABC for a general RxC rectangular matrix class.
 * The elements of the matrix are stored in a one-dimensional double
 * array in row major order. A second array of double pointers, 
 * Rows, holds pointers to each row. The element array itself is pointed to
 * by Rows[0]. Derived from Matbase_.
 */
class Rectbase_ : public virtual Matbase_
{
    // methods interface
    public:
    
	// constructors
    /* Allocates a Row x Col matrix. Without any arguments, a 3x3
     * matrix will be created. With one argument, the Col parameter
     * defaults to 0,  which is interpreted as an instruction to create
     * a square Row x Row matrix. Square matrices may be created
     * by Rectbase_(N, N) or Rectbase_(N) or Rectbase_(N, 0).
     * All elements will be initialised to 0.0.
     */
    Rectbase_(unsigned int Row=3, unsigned int Col=0);
    
    /* Allocate a Row*Col matrix as above and initialise it with
     * a traditional double ** style array Arr. It is the caller's 
     * responsibility to make sure that Arr has the correct dimensions.
     * If Arr==NULL then the matrix elements will be set to 0.0
     */
    Rectbase_(const double **Arr, unsigned int Row, unsigned int Col);
    
    /* the copy constructor */
    Rectbase_(const Rectbase_& Mat);
    
	// destructor
    virtual ~Rectbase_() { delete [] Rows[0]; }
    
	// size
    /* max_size(): returns the larger of the row or col. number. */
    unsigned int max_size() const { return((rno()>cno())? rno(): cno()); }
    
    	// access
    /* Safe access: the overlaid function call operator is used for this.
     * Error messages are printed if one of the two indices is out of
     * range and the index in question will be replaced by 0.
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
    
    /* get_array(): converts the calling object to a traditional
     * double ** array and sets Row and Col to its corresponding
     * dimensions. Returns a NULL ptr if allocation failed.
     */
    double** get_array(unsigned int& Row, unsigned int&Col) const;
        
	// end of public interface
    protected:

	// printing
    void print_rows(ostream& Out, unsigned int Sizew, unsigned int Jbeg, 
	unsigned int Items, unsigned int Width, unsigned int Prec) const;
	
	// memory management
    void init_rowptrs(unsigned int Colno);
    
};
// END OF CLASS Rectbase_

// ==== END OF HEADER Rectbase.h ====

#endif
