#ifndef FAKEBETA_CLASS
#define FAKEBETA_CLASS

// ==== PROJECT DRAGON:	HEADER Fakebeta.h ====

/* Distances between the C-alphas on the backbone and the
 * fake C-beta atoms (representing the side-chains).
 * The C-beta positions are determined by the backbone.
 */

// SGI C++ 4.0, IRIX 5.3, 28. Nov. 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Sqmat.h"
#include "Trimat.h"
#include "Points.h"

// ---- INCLUDE FILES ----

#include "Polymer.h"

// ==== CLASSES ====

/* Class Fakebeta_: stores two matrices, Distab and Distb. The former
 * holds the squared distances between the i-th CA and the j-th
 * CB; the latter holds the CB:CB squared distances. The CA:CA distance matrix
 * is not stored here.
 * The basic geometry is explained below (comment coming from
 * DRAGON 3.x):-
 * The chain is made up of a C-alpha backbone and is decorated
 * by fake C-beta atoms which represent the centroids of the
 * side chains. The monomers are 2D and therefore achiral in >=3D.
 * The i-th fake beta atom sits on the line connecting
 * the i-th C-alpha and the midpoint between the (i-1)th and
 * (i+1)th C-alphas. The first [0] and last [Rno+1] points on the backbone
 * correspond to the terminal NH3+ and COO-, respectively.
 * When I derived the distance equations, I used letters to
 * mark the atoms instead of indexing to make life easier.
 * Here is the scheme:
 * 
 *       J <-- the i-th fake C-beta atom
 *       |
 *       | <-- Dbj, the prescribed alpha:beta distance (0 for Gly)
 *       |
 *       B <-- the i-th C-alpha atom
 *     / : \
 *    /  :  \
 *   A...H...C <-- the (i+1)th C-alpha atom
 *   !   !
 *   !   ---------- the midpoint between A and C
 *   -------------- the (i-1)th C-alpha atom
 * 
 * Of course AC is not orthogonal to BH if AB != BC,  but BHJ
 * are always collinear. If ABC are collinear, then BH==0, 
 * and then J is assumed to be equivalent to B (the C-beta
 * "riding" on the C-alpha). 
 * The quantity Lambda is defined as BJ/HJ,  0...1, and B
 * divides the JH segment as BJ:BH=Lambda:(1-Lambda).
 * Lambdas can be const accessed (for the benefit of the steric routines).
 */
class Fakebeta_
{
    // data
    private:
    
    Sqmat_ Distab;  // [i][j]==dist^2(CA(i),CB(j))
    Trimat_ Distb;  // CB:CB distances
    Array_<float> Lambda, Dhj; // auxiliary vectors (see comments above)
    
    // methods
    public:
    
	// constructors
    /* Init to hold Rno monomers plus the two terminals. */
    Fakebeta_(unsigned int Rno=1): 
	Distab(Rno+2), Distb(Rno+2), Lambda(Rno+2), Dhj(Rno+2) {}
    
	// access
    double ab(unsigned int i, unsigned int j) const { return(Distab(i, j)); }
    double bb(unsigned int i, unsigned int j) const { return(Distb(i, j)); }
    const float& lambda(unsigned int i) const { return(Lambda[i]); }
    
	// distance update
    /* update(): updates the CA:CB and CB:CB distance matrices using the
     * CA:CA matrix Dista and the prescribed CA(i):CB(i) distances from
     * Polymer_ . The matrices within may be resized if necessary.
     * Return value: the new size.
     */
    unsigned int update(const Trimat_& Dista, const Polymer_& Polymer);
    
    /* beta_xyz(): generates the fake C-beta coordinates from the C-alpha coordinates
     * stored in Xyz and puts the result into Beta.
     */
    static void beta_xyz(const Points_& Xyz, const Polymer_& P, Points_& Beta);
    
    // private methods
    private:
    
    static float get_dist(float D1, float D2, float D3, float L);
    void make_lambda(const Trimat_& Dista, const Polymer_& Polymer);
};
// END OF CLASS Fakebeta_

// ==== END OF HEADER Fakebeta.h ====
#endif	/* FAKEBETA_CLASS */
