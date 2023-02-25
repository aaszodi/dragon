#ifndef HIROT_CLASS
#define HIROT_CLASS

// ==== HEADER Hirot.h ====

/* Contains an implementation of McLachlan's RMS rotation algorithm
 * (aka "Procrustes rotation" in statistics)
 * adapted to high-dimensional problems. 
 */

// SGI C++ 4.0, IRIX 5.3, 27. Oct. 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Sqmat.h"
#include "Points.h"
#include "Svd.h"

// ==== CLASSES ====

/* Hirot_: a class for performing high-dimensional rotations.
 * The object is initialised appropriately whenever a rotation
 * is needed and then can be queried for the rotation matrix
 * or the RMS value.
 */
class Hirot_
{
    // data
    private:
    
    Svd_ Svd;	// the singular value decomposition object
    Sqmat_ Mixtensor, Rot;  // the "mixed tensor" and the rotation matrix
    int Rank;	// the rank of Mixtensor (-1 if no decomp was done)
    
    // methods
    public:
    
	// constructors
    /* Default: just calls the other default ctors. */
    Hirot_(): Svd(), Mixtensor(), Rot(), Rank(-1) {}
    
	// access
    /* rot_matrix(): returns a reference to the best rotation matrix. */
    const Sqmat_& rot_matrix() const { return(Rot); }
    
    /* det_sign(): returns the sign of the determinant of the Mixtensor.
     * Returns 1 for "pure" rotations, -1 for "improper" rotations
     * (i. e. with inversion), 0 if the structures are "flat" and
     * the tensor is rank deficient. Conditions the SVD of Mixtensor
     */
    int det_sign();
    
	// rotations
    /* best_rot(): The RMS rotation algorithm adapted to
     * >=3D spaces. The point set X is rotated onto Y. Both are assumed to have
     * the same dimensionality and the same number of active points.
     * The activity masks need not be equal, but the first active
     * point in X will correspond to the first active in Y etc.
     * The active points are assumed to have been centered before the call.
     * W is a weight vector (assumed to be all-equal if omitted).
     * NOTE:  <3D point sets are not treated correctly.
     * Return value: the sign of the determinant of the transform matrix.
     */
    int best_rot(const Points_& X, const Points_& Y, const Vector_& W);
    int best_rot(const Points_& X, const Points_& Y);
    
    /* best_rotflip(): the same as hi_rot() but there's no check for
     * the det sign. This routine constructs a unitary transform matrix that
     * flips the conformations if necessary to achieve optimal superposition.
     */
    void best_rotflip(const Points_& X, const Points_& Y, const Vector_& W);
    void best_rotflip(const Points_& X, const Points_& Y);
    
	// RMS and transformation
    /* get_rms(): returns the weighted RMS value between the point sets
     * Y and Rot*X (Rot is in the calling object). 
     * Both sets should have the same no. of active points
     * and dimensions. W holds the weights for the pairs (if not uniform).
     * Returns -1.0 on error.
     */
    double get_rms(const Points_& X, const Points_& Y, const Vector_& W) const;
    double get_rms(const Points_& X, const Points_& Y) const;
    
    /* apply_transform(): replaces X with Rot*X and returns the
     * weighted RMS value between Y and Rot*X or -1.0 on error.
     * If W is not specified then uniform weights are used.
     * Rot is inside the calling object.
     */
    double apply_transform(Points_& X, const Points_& Y, const Vector_& W) const;
    double apply_transform(Points_& X, const Points_& Y) const;
    
    // hidden methods
    private:
    
    int get_rot(unsigned int Dim);
    void get_rotflip(unsigned int Dim);
    
    static int check_data(const Points_& X, const Points_& Y);
    static int check_data(const Points_& X, const Points_& Y, const Vector_& W);
    
    void make_mixtensor(const Points_& X, const Points_& Y);
    void make_mixtensor(const Points_& X, const Points_& Y, const Vector_& W);
    
};
// END OF CLASS Hirot_

// ==== END OF HEADER Hirot.h ====

#endif	/* HIROT_CLASS */

