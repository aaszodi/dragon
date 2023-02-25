#ifndef POINTS_CLASS
#define POINTS_CLASS

// ==== HEADER Points.h ====

/* A maskable array of Vector_ objects for storing and
 * manipulating point coordinates. Essentially,  an improved
 * version of the Pset_ class.
 */

// SGI C++, 30. Apr. 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>

// ---- INCLUDE FILES ----

#include "Maskarr.h"	// ancestor: maskable array template class
#include "Vector.h"	// the Vector_ class
#include "Sqmat.h"	// square matrix class
#include "Trimat.h"	// triangular matrix class

// ==== CLASSES ====

/* Class Points_ : derived from the maskable array template class
 * Maskarr_, this class stores Vector_ objects (double-precision
 * vectors) in a maskable array. Items can be activated and
 * deactivated, and access (via indexing []) is permitted to
 * active items only. The vectors may have different dimensions
 * but the active vectors' dimensions may be set to a common value
 * which is necessary for some of the "overall" operations.
 */
class Points_: public Maskarr_<Vector_>
{
    // methods
    public:
    
	// constructors and dtor
    /* Initialises the array to hold N (default=1) vectors, each
     * D-dimensional (default=3). All vectors will be active.
     */
    Points_(unsigned int N=1, unsigned int D=3);
	
    /* Initialises the array with the bitmap in Initmask. It will have 
     * Initmask.len() items having D (default=3) dimension each.
     * Activation pattern is as in Initmask, of course.
     */
    Points_(const Bits_& Initmask, unsigned int D=3);
    
    virtual ~Points_() {}
    
	// dimension
    /* dim_range(Low,High): returns the smallest and largest dimension
     * of the active vectors. Both are set to 0 if there are no active 
     * items.
     */
    void dim_range(unsigned int& Low, unsigned int& High) const;
	
    /* dim_low(), dim_high(): return the lowest and highest dimension 
     * of the set of active points or 0 if ther were no active points.
     * (If both are needed, the dim_range() method above would be slightly
     * more efficient.)
     */
    unsigned int dim_low() const;
    unsigned int dim_high() const;

    /* dim(): returns the dimension of the active vectors if they have the
     * same dimension or 0 if the dimensions are different or there are no
     * active vectors.
     * dim(D): sets the dimension of all active vectors to D. Returns old dimension
     * or 0 if they were different or there are no active vectors.
     */
    unsigned int dim() const;
    unsigned int dim(unsigned int D);
    
    /* len_dim(): adjusts the overall length to L, switches on all vectors
     * and sets their dimensions to D (default 3).
     */
    void len_dim(unsigned int L, unsigned int D=3) { len(L); mask(true); dim(D); }
    
	// arithmetics
    /* Points*=Scalar: multiplies all active vectors by Scalar in place.
     * Returns calling object.
     */
    Points_& operator*=(double Scalar);
    
    /* Points*=Matrix: premultiplies all active vectors by a square matrix
     * Matrix in place. Warnings are printed in case of dimension mismatches
     * and only those vectors which have matching dimensions will be modified.
     * Returns calling object.
     */
    Points_& operator*=(const Sqmat_& Matrix);

    /* Points+=Vector, Points-=Vector: adds/subtracts Vector to all active
     * vectors. If the dimensions don't match, then nothing happens (some
     * warnings are printed). Returns calling object.
     */
    Points_& operator+=(const Vector_& Vector);
    Points_& operator-=(const Vector_& Vector);
    
    /* centroid(W): calculates the weighted centroid of the active points,
     * W should be at least as long as the number of currently active points.
     * centroid(): If the argument is omitted, then uniform weighting is assumed.  *  * Since the dimensions need not be equal, first the maximal dimension is obtained
     * and then all active vectors in the set are temporarily "promoted"
     * to this dimension by padding them with 0-s at the end (similar to the
     * Vector_::dim(D) method). If there were no active points then a 
     * warning is printed and a 3D null-vector returned. Otherwise the 
     * centroid is returned (with the maximal dimension).
     */
    Vector_ centroid(const Vector_& W) const;
    Vector_ centroid() const;
    
    // distance matrices
    
    /* dist_mat(Dist), dist_mat2(Dist2): construct the interpoint distance
     * and squared distance matrices using the active points only. The matrices
     * will not be touched if there are no active points or if the active
     * point dimensions don't match. The matrix size will be adjusted
     * silently if necessary. Return calling object.
     */
    const Points_& dist_mat(Trimat_& Dist) const;
    const Points_& dist_mat2(Trimat_& Dist2) const;
    
    // output
    /* pdb_list(): produces a very crude PDB listing from Points provided it is
     * 3-dimensional. Amino acid code is 'G', one CA-atom per vector.
     * Return value: 0 if the dimension is outside the range [1..3], 
     * the actual dimension otherwise.
     */
    int pdb_list(ostream& Out) const;

};
// END OF CLASS Points_

// ---- PROTOTYPES ----

ostream& operator<<(ostream& Out, const Points_& Points);

// ==== END OF HEADER Points.h ====

#endif	/* POINTS_CLASS */
