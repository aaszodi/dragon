#ifndef __PSET_CLASS__
#define __PSET_CLASS__

// ==== HEADER Pset.h ====

/* The Point Set Class. A point set is an array of vectors
 * representing points in an Euclidean space. The points can be
 * switched on/off. Useful for storing atomic coordinates.
 */

// SGI C++ 3.2.1, IRIX 5.2, 21. Apr. 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- INCLUDE FILES ----

#include "Array.h"
#include "Vector.h"
#include "Bits.h"
#include "Matrix.h"
#include "Sqmat.h"
#include "Trimat.h"

// ==== CLASSES ====

/* Class Pset_: implements a point set class as an array of
 * vectors representing points in a Dim-dimensional space.
 * Each point can be "active" or "passive" i.e. switched 
 * on or off. The activity state is stored in a bit vector
 * which has the same length as the point array. The vectors
 * can be scaled (multiplied by a scalar), translated (a vector
 * can be added or subtracted) and pre-multiplied by a square matrix.
 */
class Pset_
{
    // data
    Array_<Vector_> Points;	// the point vectors
    unsigned int Dim;	    // all must have the same dimension
    Bits_ Active;   // active/passive flags
    
    public:
    // methods
    
	// constructors
    /* Init to contain N points (default 1) of D dimensions (default 3).
     * All points will be active.
     */
    Pset_(unsigned int N=1, unsigned int D=3);
    
    /* Init by an NxD matrix Mat to hold N points of D dimensions.
     * The point coordinates are assumed to be held in the rows
     * of the matrix. All points will be active.
     */
    Pset_(const Matrix_& Mat);
    
	// matrix conversion
    /* Matrix_(): converts a point set to an NxD rectangular matrix
     * with all points as rows.
     */
    operator Matrix_() const;
    
	// point access
    
    /* [] unsafe access: the index Idx is not tested. */
    const Vector_& operator[](unsigned int Idx) const { return(Points[Idx]); }
    Vector_& operator[](unsigned int Idx) { return(Points[Idx]); }
    
    /* () Safe access: the index Idx is tested. */
    const Vector_& operator()(unsigned int Idx) const { return(Points(Idx)); }
    Vector_& operator()(unsigned int Idx) { return(Points(Idx)); }
    
	// size access
    
    /* len(): returns the number of points in the object.
     * len(Size): sets the number of points to Size. If Size==0
     * then no action is taken. The new points are switched OFF
     * by default. Returns old size
     */
    unsigned int len() const { return(Points.len()); }
    unsigned int len(unsigned int Size);
    
    /* dim(): returns the current dimensionality of the points.
     * dim(D): sets the dimension of ALL (not only of the active)
     * points to D. Returns old dimension. D==0 is silently ignored.
     */
    unsigned int dim() const { return(Dim); }
    unsigned int dim(unsigned int D);

	// activation

    /* active(): returns a bit-vector with the activation flags.
     * active(Flags): sets the activation flags to Flags if the
     * dimensions match and returns the old flags.
     * active(Val): set all activation values to Val and return old flags.
     */
    Bits_ active() const { return(Active); }
    Bits_ active(const Bits_& Flags)
    {
	Bits_ Oldflags=Active;
	if (Flags.len()==Active.len()) Active=Flags;
	return(Oldflags);
    }
    Bits_ active(bool Val)
    {
	Bits_ Oldflags=Active;
	Active.set_values(Val);
	return(Oldflags);
    }
    
    /* active_no(): returns the number of active points in the set. */
    unsigned int active_no() const { return(Active.on_no()); }
    
    /* flag(Idx): returns the activation status of the Idx-th vector
     * or false if Idx is out of range.
     * flag(Idx, Val): set the activation status of the Idx-th vector
     * to Val (true/false). Do nothing for invalid Idx.
     */
    bool flag(unsigned int Idx) const { return(Active.get_bit(Idx)); }
    bool flag(unsigned int Idx, bool Val) 
    { return(Active.set_bit(Idx, Val)); }
    
	// arithmetics
    
    /* Pset*=S: scaling by the scalar S. Returns the calling object */
    Pset_& operator*=(double S);
	
    /* centroid(): calculates and returns the centroid of the active
     * points. All points have equal weights. If no active points were
     * found then the null-vector will be returned.
     */
    Vector_ centroid() const;
    
    /* Pset+=Vec: translates the active points in Pset by a vector Vec.
     * Nothing happens if Vec has wrong dimension.
     * Returns the calling object.
     */
    Pset_& operator+=(const Vector_& Vec);
    
    /* Pset-=Vec: translates the active points in Pset by a vector -Vec.
     * (Centers the point set on Vec in other words.)
     * Returns the calling object.
     */
    Pset_& operator-=(const Vector_& Vec);
    
    /* Pset*=Sqmat: premultiplies the active points in Pset by a
     * square matrix. (It would be cumbersome to define a Sqmat*Pset
     * function instead.) Nothing happens if the dimensions don't match.
     * Returns calling object.
     */
    Pset_& operator*=(const Sqmat_& Sqmat);

    /* dist_mat(),dist_mat2(): calculate ALL the interpoint distances and the
     * squared interpoint distances, respectively. Dmat is a Trimat_
     * the size of which will be adjusted silently to the number of 
     * ALL points. 
     */
    void dist_mat(Trimat_& Dmat) const;
    void dist_mat2(Trimat_& Dmat) const;
    
    friend ostream& operator<<(ostream& Out, const Pset_& Pset);
};
// END OF CLASS Pset_

// ==== END OF HEADER Pset.h ====

#endif

