#ifndef __SQBASE_CLASS__
#define __SQBASE_CLASS__

// ==== HEADER Sqbase.h ====

/* Abstract base class for square (and triangular) matrices. */

// SGI C++ 4.0, IRIX 5.3, 3. May 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

// ---- INCLUDE FILES ---- 

#include "Matbase.h"
#include "Vector.h"

// ==== CLASSES ====

/* class Sqbase_ : abstract base class for square and triangular matrices.
 * Does not add new data members. Defines main-diagonal related methods
 * and other specialities. Derived from Matbase_
 */
class Sqbase_: virtual public Matbase_
{
    // methods
    public:
    
	// empty destructor
    virtual ~Sqbase_() {}
    
	// size
    unsigned int cno() const { return(R); } // Row==Col for these
    virtual void set_size(unsigned int Size) =0;
    
	// main diagonal routines
    /* diag(): copies the diagonal into a Vector_ object.
     * diag(Vec): sets the diagonal to the values in the Vector_ object Vec.
     * If the dimensions don't match,  no action is taken. 
     */
    Vector_ diag() const;
    void diag(const Vector_& Vec);
    
    /* diag_matrix(): turns the calling object into a diagonal matrix 
     * with the value Dval in all diagonal positions. Dval==1 by default
     * so diag_matrix() produces the unit matrix w/o a parameter.
     */
    void diag_matrix(double Dval=1);
    
    /* get_trace(): returns the sum of the diagonal elements. */
    double get_trace() const;

	// transposition
    virtual void transpose_inplace() =0;
};
// END OF CLASS Sqbase_

// ==== END OF HEADER Sqbase_ ====
#endif
