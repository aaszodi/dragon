// ==== METHODS Sqbase.c++ ====

/* Abstract base class for square (and triangular) matrices. */

// SGI C++ 4.0, IRIX 5.3, 3. May 1995. Andris Aszodi

// ---- HEADER ----

#include "Sqbase.h"

// ==== METHODS ====

// ---- Diagonal access ----

/* diag(): copies the diagonal into a Vector_ object.
 * diag(Vec): sets the diagonal to the values in the Vector_ object Vec.
 * If the dimensions don't match,  no action is taken. 
 */
Vector_ Sqbase_::diag() const
{
    Vector_ Diag(rno());
    for (register unsigned int i=0; i<rno(); i++) Diag[i]=Rows[i][i];
    return(Diag);
}
// END of diag()

void Sqbase_::diag(const Vector_& Vec)
{
    if (rno()!=Vec.dim()) { prt_err(DIM_MISMATCH, "diag(V)"); return; }
    for (register unsigned int i=0; i<rno(); i++) Rows[i][i]=Vec[i];
    return;
}
// END of diag(Vec)

/* diag_matrix(): turns the calling object into a diagonal matrix 
 * with the value Dval in all diagonal positions. Dval==1 by default
 * so diag_matrix() produces the unit matrix w/o a parameter.
 */
void Sqbase_::diag_matrix(double Dval)
{
    set_values();   // zeroes
    for (register unsigned int i=0; i<rno(); i++) Rows[i][i]=Dval;
}
// END of diag_matrix()

/* get_trace(): returns the sum of the diagonal elements. */
double Sqbase_::get_trace() const
{
    double Trace=0.0;
    for (register unsigned int i=0; i<rno(); i++) Trace+=Rows[i][i];
    return(Trace);
}
// END of get_trace()

// ==== END OF METHODS Sqbase.c++ ====
